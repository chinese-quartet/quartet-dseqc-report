#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import io
import json
import logging
import os
import re
import zipfile
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff

from multiqc import config
from multiqc.plots import table, linegraph
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import report
from multiqc.modules.qualimap import QM_BamQC
from quartet_dnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Post-alignment Quality Control',
            target='postalignment_qc',
            anchor='postalignment_qc',
            href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is an report module to show the data quality after alignment.'
        )

        # Find and load any input files for postalignment_qc
        table_summary = []
        for f in self.find_log_files('postalignment_qc/postalignment_qc_summary'):
            lines = f['f'].splitlines()
            keys = lines[0].split('\t')
            content = lines[1:]
            for values in content:
                table_summary.append(dict(zip(keys, values.split('\t'))))
		
        table_summary_dic = {}
        for i in table_summary:
            key = i['Sample']
            pop_i = i.pop('Sample')
            table_summary_dic[key] = i

        if len(table_summary_dic) != 0:
            self.plot_summary_table('postalignment_qc_summary', table_summary_dic)
        else:
            log.debug('No file matched: postalignment_qc - postalignment_qc_summary.txt')	

        # Set up class objects to hold parsed data()
        self.general_stats_data = defaultdict(lambda: dict())

        # General stats - genome_results.txt
        self.qualimap_bamqc_genome_results = dict()
        for f in self.find_log_files('postalignment_qc/bamqc/genome_results'):
            QM_BamQC.parse_genome_results(self, f)
        self.qualimap_bamqc_genome_results = self.ignore_samples(self.qualimap_bamqc_genome_results)
        if len(self.qualimap_bamqc_genome_results) > 0:
            self.write_data_file(self.qualimap_bamqc_genome_results, 'multiqc_qualimap_bamqc_genome_results')

        # Coverage - coverage_histogram.txt
        self.qualimap_bamqc_coverage_hist = dict()
        for f in self.find_log_files('postalignment_qc/bamqc/coverage', filehandles=True):
            QM_BamQC.parse_coverage(self, f)
        self.qualimap_bamqc_coverage_hist = self.ignore_samples(self.qualimap_bamqc_coverage_hist)

        # Insert size - insert_size_histogram.txt
        self.qualimap_bamqc_insert_size_hist = dict()
        for f in self.find_log_files('postalignment_qc/bamqc/insert_size', filehandles=True):
            QM_BamQC.parse_insert_size(self, f)
        self.qualimap_bamqc_insert_size_hist = self.ignore_samples(self.qualimap_bamqc_insert_size_hist)

        # GC distribution - mapped_reads_gc-content_distribution.txt
        self.qualimap_bamqc_gc_content_dist = dict()
        self.qualimap_bamqc_gc_by_species = dict()  # {'HUMAN': data_dict, 'MOUSE': data_dict}
        for f in self.find_log_files('postalignment_qc/bamqc/gc_dist', filehandles=True):
            QM_BamQC.parse_gc_dist(self, f)
        self.qualimap_bamqc_gc_content_dist = self.ignore_samples(self.qualimap_bamqc_gc_content_dist)
        self.qualimap_bamqc_gc_by_species = self.ignore_samples(self.qualimap_bamqc_gc_by_species)

        num_parsed = max(
            len(self.qualimap_bamqc_genome_results),
            len(self.qualimap_bamqc_coverage_hist),
            len(self.qualimap_bamqc_insert_size_hist),
            len(self.qualimap_bamqc_gc_content_dist)
        )
        # Go no further if nothing found
        if num_parsed == 0:
            return 0

        try:
            covs = config.qualimap_config['general_stats_coverage']
            assert type(covs) == list
            assert len(covs) > 0
            covs = [str(i) for i in covs]
            log.debug('Custom Qualimap thresholds: {}'.format(', '.join([i for i in covs])))
        except (AttributeError, TypeError, AssertionError):
            covs = [1, 5, 10, 30, 50]
            covs = [str(i) for i in covs]
            log.debug('Using default Qualimap thresholds: {}'.format(', '.join([i for i in covs])))
        self.covs = covs

        # Make the plots for the report
        if len(self.qualimap_bamqc_coverage_hist)>0 and len(self.qualimap_bamqc_insert_size_hist)>0 and len(self.qualimap_bamqc_gc_content_dist)>0:
            QM_BamQC.report_sections(self)

    # Helper functions
    def get_s_name(self, f):
        s_name = os.path.basename(os.path.dirname(f['root']))
        s_name = self.clean_s_name(s_name, f['root'])
        if s_name.endswith('.qc'):
            s_name = s_name[:-3]
        return s_name

    def plot_summary_table(self, id, data, title='Table summary', section_name='Table summary', description=None, helptext=None):
        """ Create the HTML for pre-alignment qc summary """
        
        headers = OrderedDict()
        headers['% GC'] = {
            'title': '% GC',
            'description': 'Average % GC Content',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        headers['Median_insert_size'] = {
            'title': 'Ins. size',
            'description': 'Median insert size',
            'min': 0,
            'scale': 'PuOr',
            'format': '{:,.0f}'
        }

        for c in [1,5,10,30]:
            headers['% Coverage_over_{}X'.format(c)] = {
                'title': '&ge; {}X'.format(c),
                'description': 'Fraction of genome with at least {}X coverage'.format(c),
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'Blues'
            }

        for c in [20,30]:
            headers['% Q{}'.format(c)] = {
                'title': '% Q{}'.format(c),
                'description': 'Percent Q{}'.format(c),
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn-rev'
            }

        headers['Median_coverage'] = {
            'title': 'Median cov',
            'description': 'Median coverage',
            'min': 0,
            'suffix': 'X',
            'scale': 'BuPu'
        }

        headers['% Mapping'] = {
            'title': '% Aligned',
            'description': '% mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'format': '{0:.2f}'
        }

        headers['% Mismatch_rate'] = {
            'title': '% Mismatch',
            'description': 'Mismatch rate',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{0:.2f}'
        }


        table_config = {
            'namespace': 'postalignment_qc_summary',
            'id': id,
            'table_title': 'Post-alignment QC Table Summary',
            'col1_header': 'Sample',
            'no_beeswarm': True,
            'sortRows': False
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'Table summary',
            anchor = id + '_anchor',
            description = description if description else '',
            helptext = helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot = table.plot(data, headers, table_config)
        )

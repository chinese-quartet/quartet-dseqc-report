#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging
import os

from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.qualimap import QM_BamQC

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
  def __init__(self):
        
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
      return None
    
    # Initialise the parent module Class object
    super(MultiqcModule, self).__init__(
      name='Post-alignment Quality Control',
      target='Post-alignment QC',
      #anchor='post_alignment_qc',
      #href='https://github.com/chinese-quartet/quartet-dseqc-report',
      info=' is an report module to show the data quality after alignment.'
    )

    # Find and load any input files for post_alignment_qc
    table_summary = []
    for f in self.find_log_files('post_alignment_qc/summary'):
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
      self.plot_summary_table('post_alignment_qc_summary', table_summary_dic)
    else:
      log.debug('No file matched: post_alignment_qc - post_alignment.txt')	

    # Set up class objects to hold parsed data()
    self.general_stats_data = defaultdict(lambda: dict())

    # General stats - genome_results.txt
    self.qualimap_bamqc_genome_results = dict()
    for f in self.find_log_files('post_alignment_qc/bamqc/genome_results'):
      QM_BamQC.parse_genome_results(self, f)
    self.qualimap_bamqc_genome_results = self.ignore_samples(self.qualimap_bamqc_genome_results)
    if len(self.qualimap_bamqc_genome_results) > 0:
      self.write_data_file(self.qualimap_bamqc_genome_results, 'multiqc_qualimap_bamqc_genome_results')

    # Coverage - coverage_histogram.txt
    self.qualimap_bamqc_coverage_hist = dict()
    for f in self.find_log_files('post_alignment_qc/bamqc/coverage', filehandles=True):
      QM_BamQC.parse_coverage(self, f)
    self.qualimap_bamqc_coverage_hist = self.ignore_samples(self.qualimap_bamqc_coverage_hist)

    # Insert size - insert_size_histogram.txt
    self.qualimap_bamqc_insert_size_hist = dict()
    for f in self.find_log_files('post_alignment_qc/bamqc/insert_size', filehandles=True):
      QM_BamQC.parse_insert_size(self, f)
    self.qualimap_bamqc_insert_size_hist = self.ignore_samples(self.qualimap_bamqc_insert_size_hist)

    # GC distribution - mapped_reads_gc-content_distribution.txt
    self.qualimap_bamqc_gc_content_dist = dict()
    self.qualimap_bamqc_gc_by_species = dict()  # {'HUMAN': data_dict, 'MOUSE': data_dict}
    for f in self.find_log_files('post_alignment_qc/bamqc/gc_dist', filehandles=True):
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
    if num_parsed != 0:
      try:
        covs = config.qualimap_config['general_stats_coverage']
        assert type(covs) == list
        assert len(covs) > 0
        covs = [str(i) for i in covs]
        log.debug('Custom Qualimap thresholds: {}'.format(', '.join([i for i in covs])))
      except (AttributeError, TypeError, AssertionError):
        covs = [1, 5, 10, 20, 30, 50]
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

  def plot_summary_table(self, id, data, title='Summary metrics', section_name='Summary metrics', description=None, helptext=None):
    """ Create the HTML for pre-alignment qc summary """
    
    headers = OrderedDict()
    headers['%Mapping'] = {
      'title': '% Aligned',
      'description': '% mapped reads',
      'max': 100,
      'min': 0,
      'scale': 'YlGn'
    }

    headers['%Mismatch Rate'] = {
      'title': '% Mismatch',
      'description': 'Mismatch rate',
      'max': 100,
      'min': 0,
      'scale': 'RdYlGn'
    }

    headers['Mendelian Insert Size'] = {
      'title': 'Ins. size',
      'description': 'Median insert size',
      'min': 0,
      'scale': 'Purples'
    }

    for c in [20,30]:
      headers['% Q{}'.format(c)] = {
        'title': '% Q{}'.format(c),
        'description': 'Percent Q{}'.format(c),
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'RdYlGn'
      }

    headers['Mean Coverage'] = {
      'title': 'Mean cov',
      'description': 'Mean coverage',
      'min': 0,
      'suffix': 'X',
      'scale': 'Blues'
    }
    
    headers['Median Coverage'] = {
      'title': 'Median cov',
      'description': 'Median coverage',
      'min': 0,
      'suffix': 'X',
      'scale': 'Blues'
    }
     
    for c in [1,5,10,20,30,50]:
      headers['PCT_{}X'.format(c)] = {
        'title': '&ge; {}X'.format(c),
        'description': 'Fraction of genome with at least {}X coverage'.format(c),
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'YlOrRd'
      }
    
    table_config = {
      'namespace': 'post_alignment_qc_summary',
      'id': id,
      'table_title': 'Post-alignment QC Table Summary',
      'col1_header': 'Sample',
      'no_beeswarm': False,
      'sortRows': False,
      'save_file': True,
      'format': '{:,.2f}'
    }

    # Add a report section with the table
    self.add_section(
      name = section_name if section_name else 'Summary metrics',
      anchor = id + '_anchor',
      description = description if description else '',
      plot = table.plot(data, headers, table_config)
    )

#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import OrderedDict
from io import StringIO
import json
import logging
import re
import zipfile
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.express as px
import plotly.figure_factory as ff

from multiqc import config
from multiqc.plots import table, linegraph
from multiqc.modules.base_module import BaseMultiqcModule
from quartet_dnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Variant Calling Quality Control',
            target='variant_calling_qc',
            anchor='variant_calling_qc',
            href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is an report module to show the variant calling quality.'
        )

        # Find and load any input files for variant calling qc summary
        table_summary = []
        for f in self.find_log_files('variant_calling_qc/variant_calling_qc_summary'):
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
            self.plot_summary_table('variant_calling_qc_summary', table_summary_dic)
        else:
            log.debug('No file matched: variant_calling_qc - variant_calling_qc_summary.txt')
        
        # Find and load any input files for benchmark score
        for f in self.find_log_files('variant_calling_qc/benchmark_score'):
            lines = f['f'].splitlines()
            parsed_data=StringIO(f['f'])
            benchmark_score_df = pd.read_csv(parsed_data, sep='\t')
            benchmark_score_df['Batch'] = benchmark_score_df['Library'] + '_' + benchmark_score_df['Type']
            benchmark_score_df.sort_values('Benchmarking_score_all', inplace=True, ascending=False)

        if len(benchmark_score_df) != 0:
            self.plot_benchmark_score('benchmark_score', benchmark_score_df)
        else:
            log.debug('No file matched: variant_calling_qc - benchmark_score.txt')
    
    # Plot summary table
    def plot_summary_table(self, id, data, title='Table summary', section_name='Table summary', description=None, helptext=None):
        """ Create the HTML for variant calling qc summary """
        
        headers = OrderedDict()
        headers['SNV_number'] = {
            'title': 'SNV Num',
            'description': 'SNV Number',
            'scale': 'BuGn',
            'format': '{:,.f}'
        }

        headers['INDEL_number'] = {
            'title': 'INDEL Num',
            'description': 'INDEL Number',
            'scale': 'BuGn',
            'format': '{:,.f}'
        }

        headers['SNV_precision'] = {
            'title': 'SNV Precision',
            'description': 'SNV Precision',
            'suffix': '%',
            'scale': 'BuGn-rev',
            'format': '{:.2f}'
        }

        headers['SNV_recall'] = {
            'title': 'SNV Recall',
            'description': 'SNV Recall',
            'suffix': '%',
            'scale': 'RdPu',
            'format': '{:.2f}'
        }

        headers['INDEL_precision'] = {
            'title': 'INDEL Precision',
            'description': 'INDEL Precision',
            'suffix': '%',
            'scale': 'PuBu',
            'format': '{:.2f}'
        }

        headers['INDEL_recall'] = {
            'title': 'INDEL Recall',
            'description': 'INDEL Recall',
            'suffix': '%',
            'scale': 'RdGy-rev',
            'format': '{:.2f}'
        }

        table_config = {
            'namespace': 'variant_calling_qc_summary',
            'id': id,
            'table_title': 'Variant Calling QC Table Summary',
            'col1_header': 'Sample',
            'no_beeswarm': True,
            'sortRows': False
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'Table summary',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )
    
    # Plot benchmark score
    def plot_benchmark_score(self, id, bench, title='Benchmark socre', section_name='Benchmark socre', description=None, helptext=None):
        """ Creat the HTML for benchmark score """

        batch = bench.shape[0]
        indicator = bench.shape[1] - 4
        x_indicator = []
        y_batch = []
        value = []
        size = []

        col_name = list(bench.columns)
        x_label = col_name[1:indicator+1]
        y_label = list(bench['Batch'])

        for i in range(batch):
            for j in range(indicator):
                y_batch.append(i+1)
                x_indicator.append(j+1)
                value.append(bench.iat[i,j+1])
                size.append(bench.iat[i,j+1]*45)

        fig = go.Figure(data=[go.Scatter(
            x=x_indicator,
            y=y_batch,
            mode='markers',
            marker=dict(
                color=value,
                size=size,
                showscale=True),
            text=value)])
        
        fig.update_layout(template='simple_white', autosize=True, height = 800)
        x_temp = list(range(indicator))
        x_vals = [i + 1 for i in x_temp]
        fig.update_xaxes(ticktext=x_label,tickvals=x_vals,tickangle=-90)
        fig['layout']['xaxis']['side'] = 'top'
        fig['layout']['xaxis']['showgrid'] = False

        y_temp = list(range(batch))
        y_vals = [i + 1 for i in y_temp]
        fig.update_yaxes(ticktext=y_label,tickvals=y_vals)

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': '',
            'auto_margin': True
        })
        
        # Add a report section with the bubble charts
        self.add_section(
            name = section_name if section_name else 'Benchmark score',
            anchor = id + '_anchor',
            description = description if description else '',
            helptext = helptext if helptext else '''
            This section shows the performance ranking of this batch sample in all reference data sets, with a total of 4 evaluation indicators(precision, recall, reproducibility, mendelian) , compared from SNV and Indel, and finally a comprehensive benchmarking score was obtained.
            ''',
            plot = html
        )
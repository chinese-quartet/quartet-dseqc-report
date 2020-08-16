#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
from io import StringIO
import numpy as np
import json
import io
import os
import tempfile
import plotly.express as px
import plotly.figure_factory as ff

from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule
from quartet_dnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Performance Assessment',
            target='performance_assessment',
            anchor='performance_assessment',
            href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is an report module to show the performance of SNV and INDEL based on reference datasets (v202007) and Quartet family-dependent built-in genetic truth.'
        )

        # Find and load any input files for SNV performance based on reference datasets (v202007)
        for f in self.find_log_files('performance_assessment/precision_recall_snv'):
            # Create df for scatter plot
            lines = f['f'].splitlines()
            parsed_data=StringIO(f['f'])
            precision_recall_snv_df = pd.read_csv(parsed_data, sep='\t')
            precision_recall_snv_df['Batch'] = precision_recall_snv_df['Type'] + '_' + precision_recall_snv_df['Library']
            # Create dic for table
            snv_precision, library = self.parse_table_data(precision_recall_snv_df, 'Precision')
            snv_recall, library = self.parse_table_data(precision_recall_snv_df, 'Recall')
            precision_recall_snv_table_dic = self.creat_table(snv_precision, snv_recall, ['Precision (%)', 'Recall (%)'], 'SNV', library)

        self.write_data_file(precision_recall_snv_table_dic, 'precison_recall_snv', sort_cols=False, data_format='txt')

        if len(precision_recall_snv_df) != 0:
            self.plot_snv_reference_datasets('precision_recall_snv', precision_recall_snv_df, precision_recall_snv_table_dic)
        else:
            log.debug('No file matched: performance_assessment - precision_recall_snv.txt')
        
        # Find and load any input files for INDEL performance based on reference datasets (v202007)
        for f in self.find_log_files('performance_assessment/precision_recall_indel'):
            # Create df for scatter plot
            lines = f['f'].splitlines()
            parsed_data=StringIO(f['f'])
            precision_recall_indel_df = pd.read_csv(parsed_data, sep='\t')
            precision_recall_indel_df['Batch'] = precision_recall_indel_df['Type'] + '_' + precision_recall_indel_df['Library']
            # Create dic for table
            indel_precision, library = self.parse_table_data(precision_recall_indel_df, 'Precision')
            indel_recall, library = self.parse_table_data(precision_recall_indel_df, 'Recall')
            precision_recall_indel_table_dic = self.creat_table(indel_precision, indel_recall, ['Precision (%)', 'Recall (%)'], 'INDEL', library)
        
        self.write_data_file(precision_recall_indel_table_dic, 'precison_recall_indel', sort_cols=False, data_format='txt')

        if len(precision_recall_indel_df) != 0:
            self.plot_indel_reference_datasets('precision_recall_indel', precision_recall_indel_df, precision_recall_indel_table_dic)
        else:
            log.debug('No file matched: performance_assessment - precision_recall_indel.txt')
        
        # Find and load any input files for SNV performance based on Quartet family-dependent built-in genetic truth
        for f in self.find_log_files('performance_assessment/mendelian_jaccard_index_snv'):
            lines = f['f'].splitlines()
            parsed_data=StringIO(f['f'])
            mendelian_jaccard_snv_df = pd.read_csv(parsed_data, sep='\t')
            mendelian_jaccard_snv_df['Batch'] = mendelian_jaccard_snv_df['Type'] + '_' + mendelian_jaccard_snv_df['Library']          
            # Create dic for table
            snv_reproducibility, library = self.parse_table_data(mendelian_jaccard_snv_df, 'Reproducibility')
            snv_mendelian, library = self.parse_table_data(mendelian_jaccard_snv_df, 'Mendelian')
            mendelian_jaccard_snv_table_dic = self.creat_table(snv_reproducibility, snv_mendelian, ['Reproducibility (%)', 'Mendelian (%)'], 'SNV', library)

        self.write_data_file(mendelian_jaccard_snv_table_dic, 'mendelian_jaccard_snv', sort_cols=False, data_format='txt')
        
        if len(mendelian_jaccard_snv_df) != 0:
            self.plot_snv_genetic_truth('mendelian_jaccard_index_snv', mendelian_jaccard_snv_df, mendelian_jaccard_snv_table_dic)
        else:
            log.debug('No file matched: performance_assessment - mendelian_jaccard_index_snv.txt')
        
        # Find and load any input files for INDEL performance based on Quartet family-dependent built-in genetic truth
        for f in self.find_log_files('performance_assessment/mendelian_jaccard_index_indel'):
            lines = f['f'].splitlines()
            parsed_data=StringIO(f['f'])
            mendelian_jaccard_indel_df = pd.read_csv(parsed_data, sep='\t')
            mendelian_jaccard_indel_df['Batch'] = mendelian_jaccard_indel_df['Type'] + '_' + mendelian_jaccard_indel_df['Library']
            # Create dic for table
            indel_reproducibility, library = self.parse_table_data(mendelian_jaccard_indel_df, 'Reproducibility')
            indel_mendelian, library = self.parse_table_data(mendelian_jaccard_indel_df, 'Mendelian')
            mendelian_jaccard_indel_table_dic = self.creat_table(indel_reproducibility, indel_mendelian, ['Reproducibility (%)', 'Mendelian (%)'], 'INDEL', library)

        self.write_data_file(mendelian_jaccard_indel_table_dic, 'mendelian_jaccard_indel', sort_cols=False, data_format='txt')

        if len(mendelian_jaccard_indel_df) != 0:
            self.plot_indel_genetic_truth('mendelian_jaccard_index_indel', mendelian_jaccard_indel_df, mendelian_jaccard_indel_table_dic)
        else:
            log.debug('No file matched: performance_assessment - mendelian_jaccard_index_indel.txt')
    
    # Parse table data
    def parse_table_data(self, df, indicator):
        mean = df.loc[df['Type']=='Query'][indicator].mean()
        sd = df.loc[df['Type']=='Query'][indicator].std()
        value = '%.2f ± %.2f' % (mean*100, sd*100)
        
        # Add batch column into full_df
        df2 = df['Sample'].str.split(r'_(LCL\d_)|(D\d_)|()\d_\d', expand=True)
        df2['Ref_batch'] = df2[0]
        full_df = pd.concat([df, df2['Ref_batch']], axis=1)
        
        # Get query library and batch
        query_library = list(df.loc[df['Type']=='Query']['Library'])[0]
        query_batch = list(full_df.loc[full_df['Type']=='Query']['Ref_batch'])[0]
        
        # Calculate the rank in reference datasets using same library
        df_same_library = full_df.loc[full_df['Library'] == query_library]
        same_library_rank = df_same_library.groupby(['Ref_batch']).mean().rank(ascending=False).reset_index()
        n = list(same_library_rank.loc[same_library_rank['Ref_batch'] == query_batch][indicator])[0]
        N = same_library_rank.shape[0]
        same_library_rank = '%i/%i' % (n,N)
        
        # Calculate the rank in all reference datasets
        all_rank = full_df.groupby(['Ref_batch']).mean().rank(ascending=False).reset_index()
        all_n = list(all_rank.loc[all_rank['Ref_batch'] == query_batch][indicator])[0]
        all_N = all_rank.shape[0]
        all_rank = '%i/%i' % (all_n,all_N)

        return([value, all_rank, same_library_rank], query_library)
    
    # Create table dic for html
    def creat_table(self, list1, list2, indicator, var_type, library):
        table_dic = OrderedDict()
        table_dic[var_type] = {indicator[0]: list1[0], indicator[1]: list2[0]}
        table_dic['All Rank'] = {indicator[0]: list1[1], indicator[1]: list2[1]}
        table_dic['%s Rank' % library] = {indicator[0]: list1[2], indicator[1]: list2[2]}
        return(table_dic)
    
    # Write table into file
    def write_data_file(self, data, fn, sort_cols=False, data_format=None):

        # Add relevant file extension to filename
        fn = '{}.{}'.format(fn, data_format)

        # JSON encoder class to handle lambda functions
        class MQCJSONEncoder(json.JSONEncoder):
            def default(self, obj):
                if callable(obj):
                    try:
                        return obj(1)
                    except:
                        return None
                return json.JSONEncoder.default(self, obj)

        # Save file
        with io.open (os.path.join(config.data_dir, fn), 'w', encoding='utf-8') as f:
            if data_format == 'json':
                jsonstr = json.dumps(data, indent=4, cls=MQCJSONEncoder, ensure_ascii=False)
                print( jsonstr.encode('utf-8', 'ignore').decode('utf-8'), file=f)
            elif data_format == 'yaml':
                yaml.dump(data, f, default_flow_style=False)
            else:
                # Default - tab separated output
                # Convert keys to strings
                data = {str(k):v for k, v in data.items()}
                # Get all headers
                h = ['']
                for sn in data.keys():
                    for k in data[sn].keys():
                        if type(data[sn][k]) is not dict and k not in h:
                            h.append(str(k))
                if sort_cols:
                    h = sorted(h)

                # Get the rows
                rows = [ "\t".join(h) ]
                for sn in sorted(data.keys()):
                    # Make a list starting with the sample name, then each field in order of the header cols
                    l = [str(sn)] + [ str(data[sn].get(k, '')) for k in h[1:] ]
                    rows.append( "\t".join(l) )

                body = '\n'.join(rows)

                print( body.encode('utf-8', 'ignore').decode('utf-8'), file=f)
    
    # Plot SNV based on reference datasets table and scatter plot
    def plot_snv_reference_datasets(self, id, fig_data, table_data, title='SNV performance based on reference datasets (v202007)', section_name='SNV performance based on reference datasets (v202007)', description=None, helptext=None):
        """ Create the HTML for SNV performance table based on reference datasets (v202007) """
        
        headers = OrderedDict()
        headers['Precision (%)'] = {
            'title': '% Precision',
            'description': '% Precision',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        headers['Recall (%)'] = {
            'title': '% Recall',
            'description': '% Recall',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        table_config = {
            'namespace': 'precison_recall_snv',
            'id': id,
            'table_title': 'SNV performance based on Quartet family-dependent built-in genetic truth',
            'col1_header': '',
            'no_beeswarm': True,
            'sortRows': False,
            'save_file': False
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'SNV performance based on reference datasets (v202007)',
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
            plot = table.plot(table_data, headers, table_config)
        )

        fig = px.scatter(fig_data, title='SNV performance based on reference datasets (v202007)', x='Recall', y='Precision', color='Batch', marginal_y='box', marginal_x='box', template='simple_white', hover_data={'Recall': ':.3f', 'Precision': ':.3f'}, height = 650)
        #fig.update_layout(legend=dict(orientation='h', yanchor='bottom', y= -0.2, xanchor='center', x = 0.4))

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })
    
        # Add a report section with the scatter plot
        self.add_section(
            name = '',
            anchor = id + '_anchor2',
            description = '',
            helptext = '',
            plot = html
        )

    # Plot INDEL based on reference datasets table and scatter plot
    def plot_indel_reference_datasets(self, id, fig_data, table_data, title='INDEL performance based on reference datasets (v202007)', section_name='INDEL performance based on reference datasets (v202007)', description=None, helptext=None):
        """ Create the HTML for INDEL performance table based on reference datasets (v202007) """
        
        headers = OrderedDict()
        headers['Precision (%)'] = {
            'title': '% Precision',
            'description': '% Precision',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        headers['Recall (%)'] = {
            'title': '% Recall',
            'description': '% Recall',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        table_config = {
            'namespace': 'precison_recall_indel',
            'id': id,
            'table_title': 'INDEL performance based on Quartet family-dependent built-in genetic truth',
            'col1_header': '',
            'no_beeswarm': True,
            'sortRows': False,
            'save_file': False
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'INDEL performance based on reference datasets (v202007)',
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
            plot = table.plot(table_data, headers, table_config)
        )

        fig = px.scatter(fig_data, title='INDEL performance based on reference datasets (v202007)', x='Recall', y='Precision', color='Batch', marginal_y='box', marginal_x='box', template='simple_white', hover_data={'Recall': ':.3f', 'Precision': ':.3f'}, height = 650)
        #width=800, height=800

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })
    
        # Add a report section with the scatter plot
        self.add_section(
            name = '',
            anchor = id + '_anchor2',
            description = '',
            helptext = '',
            plot = html
        )

    # Plot SNV based on Quartet family-dependent built-in genetic truth table and scatter plot
    def plot_snv_genetic_truth(self, id, fig_data, table_data, title='SNV performance based on Quartet family-dependent built-in genetic truth', section_name='SNV performance based on Quartet family-dependent built-in genetic truth', description=None, helptext=None):
        """ Create the HTML for SNV performance table based on Quartet family-dependent built-in genetic truth """

        headers = OrderedDict()
        headers['Reproducibility (%)'] = {
            'title': '% Reproducibility',
            'description': '% Reproducibility',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        headers['Mendelian (%)'] = {
            'title': '% Mendelian',
            'description': '% Mendelian',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        table_config = {
            'namespace': 'mendelian_jaccard_index_snv',
            'id': id,
            'table_title': 'SNV performance based on Quartet family-dependent built-in genetic truth',
            'col1_header': '',
            'no_beeswarm': True,
            'sortRows': False,
            'save_file': False
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'SNV performance based on Quartet family-dependent built-in genetic truth',
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
            plot = table.plot(table_data, headers, table_config)
        )

        fig = px.scatter(fig_data, title='SNV performance based on Quartet family-dependent built-in genetic truth', x='Reproducibility', y='Mendelian', color='Batch', marginal_y='box', marginal_x='box', template='simple_white', hover_data={'Reproducibility': ':.3f', 'Mendelian': ':.3f'}, height = 650)

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })
    
        # Add a report section with the scatter plot
        self.add_section(
            name = '',
            anchor = id + '_anchor2',
            description = '',
            helptext = '',
            plot = html
        )        

    # Plot INDEL based on Quartet family-dependent built-in genetic truth table and scatter plot
    def plot_indel_genetic_truth(self, id, fig_data, table_data, title='INDEL performance based on Quartet family-dependent built-in genetic truth', section_name='INDEL performance based on Quartet family-dependent built-in genetic truth', description=None, helptext=None):
        """ Create the HTML for INDEL performance table based on Quartet family-dependent built-in genetic truth """

        headers = OrderedDict()
        headers['Reproducibility (%)'] = {
            'title': '% Reproducibility',
            'description': '% Reproducibility',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        headers['Mendelian (%)'] = {
            'title': '% Mendelian',
            'description': '% Mendelian',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.2f}'
        }

        table_config = {
            'namespace': 'mendelian_jaccard_index_indel',
            'id': id,
            'table_title': 'INDEL performance based on Quartet family-dependent built-in genetic truth',
            'col1_header': '',
            'no_beeswarm': True,
            'sortRows': False,
            'save_file': False
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'INDEL performance based on Quartet family-dependent built-in genetic truth',
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
            plot = table.plot(table_data, headers, table_config)
        )

        fig = px.scatter(fig_data, title='INDEL performance based on Quartet family-dependent built-in genetic truth', x='Reproducibility', y='Mendelian', color='Batch', marginal_y='box', marginal_x='box', template='simple_white', hover_data={'Reproducibility': ':.3f', 'Mendelian': ':.3f'}, height = 650)

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the scatter plot
        self.add_section(
            name = '',
            anchor = id + '_anchor2',
            description = '',
            helptext = '',
            plot = html
        )   


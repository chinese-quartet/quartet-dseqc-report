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
                
        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None
        
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Variant Calling Quality Control',
            target='variant_calling_qc',
            anchor='variant_calling_qc',
            href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is an report module to show quality assessment of the variant calling.'
        )

        # Find and load input files
        ## historical batches' performance     
        for f in self.find_log_files('variant_calling_qc/history'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            history_df = pd.read_csv(f_p, sep='\t')
        if len(history_df) == 0:
            log.debug('No file matched: variant_calling_qc - history.txt')
        
        ## snv_indel_summary (variants.calling.qc.txt)
        for f in self.find_log_files('variant_calling_qc/snv_indel_summary'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            snv_indel_df = pd.read_csv(f_p, sep='\t')
        # Plot the variant calling qc summary table
        if len(snv_indel_df) == 0:
            log.debug('No file matched: variant_calling_qc - variants.calling.qc.txt')
        else:
            snv_indel_summary = []
            df = snv_indel_df
            # Transfer into the required format
            keys = df.columns.tolist()
            for index in df.index:
                row = df.loc[index].tolist()
                snv_indel_summary.append(dict(zip(keys, row)))
            
            snv_indel_summary_dic = {}
            for i in snv_indel_summary:
                key = i['Sample']
                pop_i = i.pop('Sample')
                snv_indel_summary_dic[key] = i
            
            self.plot_table_based_on_reference_datasets('snv_indel_summary', snv_indel_summary_dic)


        ## mendelian_summary (mendelian.txt)
        for f in self.find_log_files('variant_calling_qc/mendelian_summary'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            mendelian_df = pd.read_csv(f_p, sep='\t')
        # Plot the variant calling qc summary table
        if len(mendelian_df) == 0:
            log.debug('No file matched: variant_calling_qc - mendelian.txt')
        else:
            df = mendelian_df
            Family = pd.DataFrame(["Family 1", "Family 2", "Family 3"], columns = ["Family"])
            # Extract SNV
            snv_tmp = df[df.Family.str.contains("SNV$")]
            snv_tmp = snv_tmp.drop(["Family"], axis = 1)
            snv_tmp = snv_tmp.reset_index(drop=True)
            snv_tmp = pd.concat([Family, snv_tmp], axis = 1)
            snv_tmp.columns = ["Family", "SNV Detected Variants", "SNV Mendelian Consistent Variants", "SNV MCR Rate"]

            # Extract INDEL
            indel_tmp = df[df.Family.str.contains("INDEL$")]
            indel_tmp = indel_tmp.drop(["Family"], axis = 1)
            indel_tmp = indel_tmp.reset_index(drop=True)
            indel_tmp = pd.concat([Family, indel_tmp], axis = 1)
            indel_tmp.columns = ["Family", "INDEL Detected Variants", "INDEL Mendelian Consistent Variants", "INDEL MCR Rate"]
            
            # Combine SNV & INDEL
            mendelian_df_comb = pd.merge(snv_tmp, indel_tmp, on = ["Family"])

            # Transfer into the required format
            mendelian_summary = []
            keys = mendelian_df_comb.columns.tolist()
            for index in mendelian_df_comb.index:
                row = mendelian_df_comb.loc[index].tolist()
                mendelian_summary.append(dict(zip(keys, row)))
            
            mendelian_summary_dic = {}
            for i in mendelian_summary:
                key = i['Family']
                pop_i = i.pop('Family')
                mendelian_summary_dic[key] = i
            
            self.plot_table_based_on_quartet_builtin_truth('mendelian_summary', mendelian_summary_dic)
    
    # Plot detailed numbers of performance assessment based on reference datasets
    def plot_table_based_on_reference_datasets(self, id, data, title='Detailed numbers of performance assessment based on reference datasets', section_name='Detailed numbers of performance assessment based on reference datasets', description="Detailed numbers of performance assessment based on reference datasets.", helptext=None):
        """ Create the HTML for detailed numbers of performance assessment based on reference datasets """
        headers = OrderedDict()
        headers['SNV number'] = {
            'title': 'SNV Num',
            'description': 'SNV Total Number',
            'scale': False,
            'format': '{:.0f}'
        }
        
        headers['SNV query'] = {
            'title': 'SNV Num in Benchmark Region',
            'description': 'SNV Nubmer in Benchmark Region',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['SNV TP'] = {
            'title': 'SNV TP',
            'description': 'SNV True Positive',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['SNV FP'] = {
            'title': 'SNV FP',
            'description': 'SNV False Positive',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['SNV FN'] = {
            'title': 'SNV FN',
            'description': 'SNV False Negative',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['INDEL number'] = {
            'title': 'INDEL Num',
            'description': 'INDEL Total Number',
            'scale': False,
            'format': '{:.0f}'
        }
        
        headers['INDEL query'] = {
            'title': 'INDEL Num in Benchmark Region',
            'description': 'INDEL Nubmer in Benchmark Region',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['INDEL TP'] = {
            'title': 'INDEL TP',
            'description': 'INDEL True Positive',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['INDEL FP'] = {
            'title': 'INDEL FP',
            'description': 'INDEL False Positive',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['INDEL FN'] = {
            'title': 'INDEL FN',
            'description': 'INDEL False Negative',
            'scale': False,
            'format': '{:.0f}'
        }

        table_config = {
            'namespace': 'snv_indel_summary',
            'id': id,
            'table_title': 'Detailed numbers of performance assessment based on reference datasets',
            'col1_header': 'Sample',
            'no_beeswarm': True,
            'sortRows': False,
            'format': '{:.0f}',
            'max_table_rows': 20,
            'decimalPoint_format': ',',
            'thousandsSep_format': ","
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'Detailed numbers of performance assessment based on reference datasets',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )
    

    # Plot detailed numbers of performance assessment based on Quartet genetic built-in truth
    def plot_table_based_on_quartet_builtin_truth(self, id, data, title='Detailed numbers of performance assessment based on Quartet genetic built-in truth', section_name='Detailed numbers of performance assessment based on Quartet genetic built-in truth', description="Detailed numbers of performance assessment based on Quartet genetic built-in truth.", helptext=None):
        """ Create the HTML for detailed numbers of performance assessment based on Quartet genetic built-in truth """
        
        headers = OrderedDict()
        headers['SNV Detected Variants'] = {
            'title': 'SNV Detected Variants',
            'description': 'SNV All Detected Variants',
            'scale': False
        }
        
        headers['SNV Mendelian Consistent Variants'] = {
            'title': 'SNV Mendelian Consistent Variants',
            'description': 'SNV Mendelian Consistent Variants',
            'scale': False
        }

        headers['INDEL Detected Variants'] = {
            'title': 'INDEL Detected Variants',
            'description': 'INDEL All Detected Variants',
            'scale': False
        }
        
        headers['INDEL Mendelian Consistent Variants'] = {
            'title': 'INDEL Mendelian Consistent Variants',
            'description': 'INDEL Mendelian Consistent Variants',
            'scale': False
        }

        table_config = {
            'namespace': 'mendelian_summary',
            'id': id,
            'table_title': 'Detailed numbers of performance assessment based on Quartet genetic built-in truth',
            'col1_header': 'Family',
            'no_beeswarm': True,
            'sortRows': False,
            'format': '{:.0f}',
            'max_table_rows': 20
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else 'Detailed numbers of performance assessment based on Quartet genetic built-in truth',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )
    

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
    def plot_snv_reference_datasets(self, id, fig_data, table_data, title='SNV performance based on reference datasets', section_name='SNV performance based on reference datasets', description=None, helptext=None):
        """ Create the HTML for SNV performance table based on reference datasets """
        
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
            name = section_name if section_name else 'SNV performance based on reference datasets',
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

        fig = px.scatter(fig_data, title='SNV performance based on reference datasets', x='Recall', y='Precision', color='Batch', template='simple_white', hover_data={'Recall': ':.3f', 'Precision': ':.3f'}, height = 650)
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
    def plot_indel_reference_datasets(self, id, fig_data, table_data, title='INDEL performance based on reference datasets', section_name='INDEL performance based on reference datasets', description=None, helptext=None):
        """ Create the HTML for INDEL performance table based on reference datasets """
        
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
            name = section_name if section_name else 'INDEL performance based on reference datasets',
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

        fig = px.scatter(fig_data, title='INDEL performance based on reference datasets', x='Recall', y='Precision', color='Batch', template='simple_white', hover_data={'Recall': ':.3f', 'Precision': ':.3f'}, height = 650)
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

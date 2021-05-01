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
from multiqc.plots import table, linegraph, scatter
from multiqc.modules.base_module import BaseMultiqcModule

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
            target='Variant calling QC',
            #anchor='Variant Calling QC',
            #href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is an report module to show quality assessment of the variant calling.'
        )

        # Find and load input files
        ## historical batches' performance     
        for f in self.find_log_files('variant_calling_qc/history'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            history_df = pd.read_csv(f_p, sep='\t')
        if len(history_df) == 0:
            log.debug('No file matched: variant_calling_qc - history.txt')
        
        ## precision_recall (reference_datasets_aver-std.txt)
        for f in self.find_log_files('variant_calling_qc/precision_recall'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            precision_recall_df = pd.read_csv(f_p, sep=' ', index_col = None, names = ["mean", "sd", "type", "indicator"])
        if len(precision_recall_df) == 0:
            log.debug('No file matched: variant_calling_qc - reference_datasets_aver-std.txt')
        else:
            df = precision_recall_df
            df["type"] = pd.DataFrame(["SNV", "INDEL", "SNV", "INDEL", "SNV", "INDEL"])
            df["indicator"] = pd.DataFrame(["Precision", "Precision", "Recall", "Recall", "F1-score", "F1-score"])

            # Get the rank of this batch
            snv_rank = self.get_snv_indel_rank(history_df, df, "SNV")
            indel_rank = self.get_snv_indel_rank(history_df, df, "INDEL")

            df["mean"] = df["mean"].round(2)
            df["sd"] = df["sd"].round(2)
            df["result"] = df["mean"].astype("string").str.cat(df["sd"].astype("string"), sep = " ± ")

            df = df[["type", "result", "indicator"]]
            snv_df = df[df["type"] == "SNV"].T
            snv_df.columns = snv_df.loc["indicator"]
            snv_df = snv_df.drop(["type", "indicator"], axis = 0)
            snv_df["Type"] = "SNV"
            snv_df["Rank"] = snv_rank

            indel_df = df[df["type"] == "INDEL"].T
            indel_df.columns = indel_df.loc["indicator"]
            indel_df = indel_df.drop(["type", "indicator"], axis = 0)
            indel_df["Type"] = "INDEL"
            indel_df["Rank"] = indel_rank

            df = pd.concat([snv_df, indel_df], axis = 0).reset_index(drop = True)

            # Transfer into the required format
            precision_recall_summary_dic = self.convert_input_data_format(df, 'Type')

            # Plot detailed numbers of performance assessment based on reference datasets
            self.assessment_based_on_reference_datasets('precision_recall_summary', precision_recall_summary_dic)

        ## quartet_mendelian (quartet_indel.txt; quartet_snv.txt)
        for f in self.find_log_files('variant_calling_qc/quartet_snv'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            quartet_snv_df = pd.read_csv(f_p, sep=' ', index_col = None, names = ["mean", "sd", "Type"])
        if len(quartet_snv_df) == 0:
            log.debug('No file matched: variant_calling_qc - quartet_snv_aver-std.txt')
        else:
            quartet_snv_df["Type"] = "SNV"
            quartet_snv_df["Rank"] = self.get_mendelian_rank(history_df, quartet_snv_df, "SNV")
            quartet_snv_df["mean"] = quartet_snv_df["mean"].round(2)
            quartet_snv_df["sd"] = quartet_snv_df["sd"].round(2)
            quartet_snv_df["MCR"] = quartet_snv_df["mean"].astype("string").str.cat(quartet_snv_df["sd"].astype("string"), sep = " ± ")

        for f in self.find_log_files('variant_calling_qc/quartet_indel'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            quartet_indel_df = pd.read_csv(f_p, sep=' ', index_col = None, names = ["mean", "sd", "Type"])
        if len(quartet_indel_df) == 0:
            log.debug('No file matched: variant_calling_qc - quartet_indel_aver-std.txt')
        else:
            quartet_indel_df["Type"] = "INDEL"
            quartet_indel_df["Rank"] = self.get_mendelian_rank(history_df, quartet_indel_df, "INDEL")
            quartet_indel_df["mean"] = quartet_indel_df["mean"].round(2)
            quartet_indel_df["sd"] = quartet_indel_df["sd"].round(2)
            quartet_indel_df["MCR"] = quartet_indel_df["mean"].astype("string").str.cat(quartet_indel_df["sd"].astype("string"), sep = " ± ")
        
        df = pd.concat([quartet_snv_df, quartet_indel_df], axis = 0).reset_index(drop = True)

        # Plot detailed numbers of performance assessment based on Quartet genetic built-in truth
        # Transfer into the required format
        quartet_mendelian_summary_dic = self.convert_input_data_format(df, 'Type')
        
        self.assessment_based_on_quartet_builtin_truth('quartet_mendelian_summary', quartet_mendelian_summary_dic)

        ## snv_indel_summary (variants.calling.qc.txt)
        for f in self.find_log_files('variant_calling_qc/snv_indel_summary'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            snv_indel_df = pd.read_csv(f_p, sep='\t')
        # Plot detailed numbers of performance assessment based on reference datasets
        if len(snv_indel_df) == 0:
            log.debug('No file matched: variant_calling_qc - variants.calling.qc.txt')
        else:
            snv_indel_summary = []
            df = snv_indel_df

            # Transfer into the required format
            snv_indel_summary_dic = self.convert_input_data_format(df, 'Sample')
        

        ## mendelian_summary (mendelian.txt)
        for f in self.find_log_files('variant_calling_qc/mendelian_summary'):
            f_p = '%s/%s' % (f['root'], f['fn'])
            mendelian_df = pd.read_csv(f_p, sep='\t')
        # Plot detailed numbers of performance assessment based on Quartet genetic built-in truth
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
            mendelian_summary_dic = self.convert_input_data_format(mendelian_df_comb, 'Family')
            

        # Figure data
        # F1-score of this batch
        f1_snv = snv_indel_df[["Sample", "SNV F1"]]
        f1_snv.columns = ["Sample", "F1-score"]
        f1_snv["Type"] = "SNV"
        snv_mcr = []
        for index, row in f1_snv.iterrows():
            if "_1_" in row.Sample:
                MCR = snv_tmp[snv_tmp.Family == "Family 1"]["SNV MCR Rate"].to_list()[0]
            elif "_2_" in row.Sample:
                MCR = snv_tmp[snv_tmp.Family == "Family 2"]["SNV MCR Rate"].to_list()[0]
            elif "_3_" in row.Sample:
                MCR = snv_tmp[snv_tmp.Family == "Family 3"]["SNV MCR Rate"].to_list()[0]
            snv_mcr.append(MCR)

        f1_snv_mcr = pd.concat([f1_snv, pd.DataFrame(snv_mcr, columns = ["MCR"])], axis = 1)
        
        f1_indel = snv_indel_df[["Sample", "INDEL F1"]]
        f1_indel.columns = ["Sample", "F1-score"]
        f1_indel["Type"] = "INDEL"
        indel_mcr = []
        for index, row in f1_indel.iterrows():
            if "_1_" in row.Sample:
                MCR = indel_tmp[indel_tmp.Family == "Family 1"]["INDEL MCR Rate"].to_list()[0]
            elif "_2_" in row.Sample:
                MCR = indel_tmp[indel_tmp.Family == "Family 2"]["INDEL MCR Rate"].to_list()[0]
            elif "_3_" in row.Sample:
                MCR = indel_tmp[indel_tmp.Family == "Family 3"]["INDEL MCR Rate"].to_list()[0]
            indel_mcr.append(MCR)

        f1_indel_mcr = pd.concat([f1_indel, pd.DataFrame(indel_mcr, columns = ["MCR"])], axis = 1)
        
        this_batch = pd.concat([f1_snv_mcr, f1_indel_mcr], axis = 0)
        this_batch["Batch"] = "Your Datasets"

        # Prepare data for scatter plot
        history_tmp_df = history_df[["sample", "mendelian", "f1", "type"]]
        history_tmp_df.columns = ["Sample", "MCR", "F1-score", "Type"]
        history_tmp_df["Batch"] = "Rest Submmited Datasets"
        # Add the F1-score and MCR info of this batch
        figure_data = pd.concat([history_tmp_df, this_batch], axis = 0)

        # SNV Performance
        snv_fig_data = figure_data[figure_data.Type == "SNV"]
        self.quartet_scatter_plot(snv_fig_data, "snv_performance", "SNV", "SNV Performance", "snv_performance")        
        
        # INDEL Performance
        indel_fig_data = figure_data[figure_data.Type == "INDEL"]
        self.quartet_scatter_plot(indel_fig_data, "indel_performance", "INDEL", "INDEL Performance", "INDEL_performance")
        
        # Detailed table 1
        self.detailed_based_on_reference_datasets('snv_indel_summary', snv_indel_summary_dic)
        # Detailed table 2
        self.detailed_based_on_quartet_builtin_truth('mendelian_summary', mendelian_summary_dic)

    # Functions for getting the rank of submitted data
    def get_snv_indel_rank(self, history_df, df, typ):
        history = history_df[history_df.type == typ].f1.to_list()
        history.sort()
        query = df[(df.indicator == "F1-score") & (df.type == typ)]["mean"].to_list()[0]
        rank = len(history) + 1
        for i in history:
            if query > i:
                rank = rank - 1
        
        rank_in_history = "%s/%s" % (rank, len(history) + 1)
        return(rank_in_history)

    def get_mendelian_rank(self, history_df, df, typ):
        history = history_df[history_df.type == typ].mendelian.to_list()
        history.sort()
        query = df[df.Type == typ]["mean"].to_list()[0]
        rank = len(history) + 1
        for i in history:
            if query > i:
                rank = rank - 1
        
        rank_in_history = "%s/%s" % (rank, len(history) + 1)
        return(rank_in_history)
    
    def convert_input_data_format(self, df, col_name):
        convert_df = []
        keys = df.columns.tolist()
        for index in df.index:
            row = df.loc[index].tolist()
            convert_df.append(dict(zip(keys, row)))
        
        convert_dic = {}
        for i in convert_df:
            key = i[col_name]
            pop_i = i.pop(col_name)
            convert_dic[key] = i
        
        return(convert_dic)

    # Functions for tables and scatter plots
    ## Assessment based on reference datasets (v202103)
    def assessment_based_on_reference_datasets(self, id, data, title='Assessment based on reference datasets', section_name='Assessment based on reference datasets', description="The reference dataset version is v202103.", helptext=None):
        """ Create the HTML for assessment based on reference datasets """
        
        headers = OrderedDict()
        headers['Precision'] = {
            'title': 'Precision',
            'description': 'Precision',
            'scale': False,
            'format': '{:.0f}'
        }
        
        headers['Recall'] = {
            'title': 'Recall',
            'description': 'Recall',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['F1-score'] = {
            'title': 'F1-score',
            'description': 'F1-score',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['Rank'] = {
            'title': 'Rank',
            'description': 'Rank',
            'scale': False,
            'format': '{:.0f}'
        }

        table_config = {
            'namespace': 'assessment_based_on_reference_datasets',
            'id': id,
            'table_title': 'Assessment based on reference datasets',
            'col1_header': '',
            'no_beeswarm': True,
            'sortRows': False,
            'format': '{:.0f}',
            'max_table_rows': 20,
            'decimalPoint_format': ',',
            'thousandsSep_format': ","
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else '',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )
    
    ## Assessment based on Quartet family-dependent built-in genetic truth
    def assessment_based_on_quartet_builtin_truth(self, id, data, title='Assessment based on Quartet family-dependent built-in genetic truth', section_name='Assessment based on Quartet family-dependent built-in genetic truth', description="", helptext=None):
        """ Create the HTML for assessment based on reference datasets """
        
        headers = OrderedDict()
        headers['MCR'] = {
            'title': 'Mendelian Concordance Rate',
            'description': 'Mendelian Consistent Rate',
            'scale': False,
            'format': '{:.0f}'
        }

        headers['Rank'] = {
            'title': 'Rank',
            'description': 'Rank',
            'scale': False,
            'format': '{:.0f}'
        }

        table_config = {
            'namespace': 'assessment_based_on_reference_datasets',
            'id': id,
            'table_title': 'Assessment based on Quartet family-dependent built-in genetic truth',
            'col1_header': '',
            'no_beeswarm': True,
            'sortRows': False,
            'format': '{:.0f}',
            'max_table_rows': 20,
            'decimalPoint_format': ',',
            'thousandsSep_format': ","
        }

        # Add a report section with the table
        self.add_section(
            name = section_name if section_name else '',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )
    
    ## Scatter plot
    def quartet_scatter_plot(self, figure_data, pconfig_id, pconfig_title, name, anchor):
        data = dict()
        for index, row in figure_data.iterrows():
            s_name = row["Sample"]
            data[s_name] = {
                'x': row["F1-score"],
                'y': row["MCR"]
            }

            if row["Batch"] == "Your Datasets":
                # blue
                data[s_name]['color'] = 'rgba(109, 164, 202, 0.9)'
            else:
                # yellow
                data[s_name]['color'] = 'rgba(250, 160, 81, 0.8)'
                # green: rgba(43, 159, 43, 0.8)
        
        pconfig = {
            'id': pconfig_id,
            'title': pconfig_title,
            'xlab': 'F1-score',
            'ylab': 'Mendelian Concordance Rate',
        }

        if len(data) > 0:
            self.add_section (
                name = name,
                anchor = anchor,
                description = """Points are coloured as follows:
                <span style="color: #6DA4CA;">Your Datasets</span>,
                <span style="color: #FAA051;">Rest Submmited Datasets</span>.""",
                plot = scatter.plot(data, pconfig)
            )

    ## Plot detailed numbers of performance assessment based on reference datasets
    def detailed_based_on_reference_datasets(self, id, data, title='Detailed numbers of performance assessment based on reference datasets', section_name='Detailed numbers of performance assessment based on reference datasets', description="", helptext=None):
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
            name = section_name if section_name else '',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )
    

    ## Plot detailed numbers of performance assessment based on Quartet genetic built-in truth
    def detailed_based_on_quartet_builtin_truth(self, id, data, title='Detailed numbers of performance assessment based on Quartet genetic built-in truth', section_name='Detailed numbers of performance assessment based on Quartet genetic built-in truth', description="", helptext=None):
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
            name = section_name if section_name else '',
            anchor = id + '_anchor',
            description = description if description else '',
            plot = table.plot(data, headers, table_config)
        )

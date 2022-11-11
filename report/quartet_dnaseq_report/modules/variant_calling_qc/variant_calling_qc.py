#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import pandas as pd
import numpy as np
import seaborn as sns

from multiqc import config
from multiqc.plots import table, scatter
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
      info=' is an report module to show quality assessment of the variant calling.'
    )
    
    # Find and load input files
    ### SUMMARY TABLE 1
    snv_indel_df = pd.DataFrame()
    for f in self.find_log_files('conclusion/precision_recall_summary'):
      if f is None:
        log.debug('No file matched: variant_calling_qc - variants.calling.qc.txt')
      else:
        f_p = '%s/%s' % (f['root'], f['fn'])
        tmp_df = pd.read_csv(f_p, sep='\t')
        tmp_df[['SNV precision', 'INDEL precision', 'SNV recall', 'INDEL recall']] = round(tmp_df[['SNV precision', 'INDEL precision', 'SNV recall', 'INDEL recall']]/100, 4)
        snv_indel_df = pd.concat([snv_indel_df, tmp_df], axis=0)
    
    snv_indel_df = snv_indel_df.reset_index(drop=True)
    # Generate F1-score
    if 'SNV F1' not in snv_indel_df.columns.to_list():
      snv_indel_df['SNV F1'] = 2*snv_indel_df['SNV precision']*snv_indel_df['SNV recall']/(snv_indel_df['SNV precision']+snv_indel_df['SNV recall'])
      snv_indel_df['INDEL F1'] = 2*snv_indel_df['INDEL precision']*snv_indel_df['INDEL recall']/(snv_indel_df['INDEL precision']+snv_indel_df['INDEL recall'])
    else:
      snv_indel_df['SNV F1'] = snv_indel_df['SNV F1']/100
      snv_indel_df['INDEL F1'] = snv_indel_df['INDEL F1']/100
    
    ### SUMMARY TABLE 2
    mendelian_df = pd.DataFrame()
    n = 1
    for f in self.find_log_files('conclusion/mendelian_summary'):
      if f is None:
        log.debug('No file matched: variant_calling_qc - project.summary.txt')
      else:
        f_p = os.path.join(f['root'], f['fn'])
        tmp_df = pd.read_csv(f_p, sep='\t')
        tmp_df.Family = 'Family %i.' % n + tmp_df.Family
        mendelian_df = pd.concat([mendelian_df, tmp_df], axis=0)
        n = n+1
   
    mendelian_df = mendelian_df.reset_index(drop=True)
    if len(mendelian_df) != 0:
      df = mendelian_df
      # Extract SNV
      snv_tmp = df[df.Family.str.contains("SNV$")].drop(["Family"], axis = 1).reset_index(drop=True)
      snv_tmp.insert(0, "Family", ["Queried_Data_Set%s" % i for i in range(1, len(snv_tmp.index)+1)], allow_duplicates=True)
      snv_mcr = "%s ± %s" % (round(snv_tmp['Mendelian_Concordance_Rate'].mean(), 2), round(np.std(snv_tmp['Mendelian_Concordance_Rate']), 2))
      # Extract INDEL
      indel_tmp = df[df.Family.str.contains("INDEL$")].drop(["Family"], axis = 1).reset_index(drop=True)
      indel_tmp.insert(0, "Family", ["Queried_Data_Set%s" % i for i in range(1, len(indel_tmp.index)+1)], allow_duplicates=True)
      indel_mcr = "%s ± %s" % (round(indel_tmp['Mendelian_Concordance_Rate'].mean(), 2), round(np.std(snv_tmp['Mendelian_Concordance_Rate']), 2))
    
    ### DETAILED TABLE 1: snv_indel_summary (variants.calling.qc.txt)
    # Plot detailed numbers of performance assessment based on reference datasets
    # Transfer into the required format
    detail1_dic = self.convert_input_data_format(snv_indel_df, 'Sample')
    self.detail_1('variant_calling_qc_details', detail1_dic)

    detail1_path = os.path.join(config.output_dir, "variant_calling_qc_details.txt")
    snv_indel_df.to_csv(detail1_path, sep="\t", index=0)

    ### DETAILED TABLE 2: mendelian_summary (mendelian.txt)
    snv_tmp.columns = ["", "SNV Detected Variants", "SNV Mendelian Consistent Variants", "SNV MCR Rate"]
    indel_tmp.columns = ["", "INDEL Detected Variants", "INDEL Mendelian Consistent Variants", "INDEL MCR Rate"]
    # Combine SNV & INDEL
    mendelian_df_comb = pd.merge(snv_tmp, indel_tmp, on = [""])
    mendelian_df_comb = mendelian_df_comb.drop(['SNV MCR Rate', 'INDEL MCR Rate'], axis=1)
    # Transfer into the required format
    detail2_dic = self.convert_input_data_format(mendelian_df_comb, "")
    self.detail_2('mendelian_details', detail2_dic)
    
    detail2_path = os.path.join(config.output_dir, "mendelian_details.txt")
    mendelian_df_comb.to_csv(detail2_path, sep="\t", index=0)
  
  ### Functions for getting the rank of submitted data
  def convert_input_data_format(self, df, col_name):
    convert_df = []
    keys = df.columns.tolist()
    for index in df.index:
      row = df.iloc[index].tolist()
      convert_df.append(dict(zip(keys, row)))
    
    convert_dic = {}
    for i in convert_df:
      key = i[col_name]
      pop_i = i.pop(col_name)
      convert_dic[key] = i
    
    return(convert_dic)
  
  ### Function 1: Plot detailed numbers of performance assessment based on reference datasets
  def detail_1(self, id, data, title='Details based on reference datasets', section_name='Details based on reference datasets', description="", helptext=None):
    """ Create the HTML for detailed numbers of performance assessment based on reference datasets """

    headers = OrderedDict()
    headers['SNV number'] = {
      'title': 'SNV num',
      'description': 'SNV total number',
      'scale': False
    }
    
    headers['SNV query'] = {
      'title': 'SNV query',
      'description': 'The intersection loci between the capture region and the high confidence interval',
      'scale': False
    }
    
    headers['SNV TP'] = {
      'title': 'SNV TP',
      'description': 'The number of detected SNVs that are consistent with reference dataset',
      'scale': False
    }
    
    headers['SNV FP'] = {
      'title': 'SNV FP',
      'description': 'The number of detected SNVs that are not mutated in the reference dataset',
      'scale': False
    }
    
    headers['SNV FN'] = {
      'title': 'SNV FN',
      'description': 'The number of SNVs not detected but validated in reference dataset',
      'scale': False
    }
    
    headers['SNV precision'] = {
      'title': 'SNV precision',
      'description': 'SNV precision',
      "min": 0,
      "max": 1,
      "scale": "Set3",
      'format': '{:,.4f}'
    }
    
    headers['SNV recall'] = {
      'title': 'SNV recall',
      'description': 'SNV recall',
      "min": 0,
      "max": 1,
      "scale": "Set2",
      'format': '{:,.4f}'
    }
    
    headers['SNV F1'] = {
      'title': 'SNV F1',
      'description': 'SNV F1',
      "min": 0,
      "max": 1,
      "scale": "Spectral",
      'format': '{:,.4f}'
    }
    
    headers['INDEL number'] = {
      'title': 'INDEL num',
      'description': 'INDEL total number',
      'scale': False
    }
    
    headers['INDEL query'] = {
      'title': 'INDEL query',
      'description': 'The intersection loci between the capture region and the high confidence interval',
      'scale': False
    }
    
    headers['INDEL TP'] = {
      'title': 'INDEL TP',
      'description': 'The number of detected INDELs that are consistent with reference dataset',
      'scale': False
    }
    
    headers['INDEL FP'] = {
      'title': 'INDEL FP',
      'description': 'The number of detected INDELs that are not mutated in the reference dataset',
      'scale': False
    }
    
    headers['INDEL FN'] = {
      'title': 'INDEL FN',
      'description': 'The number of INDELs not detected but validated in reference dataset',
      'scale': False
    }
    
    headers['INDEL precision'] = {
      'title': 'INDEL precision',
      'description': 'INDEL precision',
      "min": 0,
      "max": 1,
      "scale": "RdYlGn-rev",
      'format': '{:,.4f}'
    }
    
    headers['INDEL recall'] = {
      'title': 'INDEL recall',
      'description': 'INDEL recall',
      "min": 0,
      "max": 1,
      "scale": "YlGnBu",
      'format': '{:,.4f}'
    }
    
    headers['INDEL F1'] = {
      'title': 'INDEL F1',
      'description': 'INDEL F1',
      "min": 0,
      "max": 1,
      "scale": "YlGn",
      'format': '{:,.4f}'
    }
    
    table_config = {
      'namespace': 'variant_calling_qc_details',
      'id': id,
      'table_title': 'Details based on reference datasets',
      'col1_header': 'Sample',
      'no_beeswarm': True,
      'sortRows': False,
      'format': '{:,.0f}',
      'max_table_rows': 20,
      'decimalPoint_format': ',',
      'thousandsSep_format': ",",
      'save_file': True
    }

    # Add a report section with the table
    self.add_section(
      name = section_name if section_name else '',
      anchor = id + '_anchor',
      description = description if description else '',
      plot = table.plot(data, headers, table_config)
    )
  

  ### Function 2: Plot detailed numbers of performance assessment based on Quartet genetic built-in truth
  def detail_2(self, id, data, title='Details based on Quartet genetic built-in truth', section_name='Details based on Quartet genetic built-in truth', description="Each row represents a set of Quartet samples, i.e. one each of D5, D6, F7 and M8. When multiple sets of technical replicates are measured, the performance of each set will be represented by row.", helptext=None):
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
      'namespace': 'mendelian_details',
      'id': id,
      'table_title': 'Details based on Quartet genetic built-in truth',
      'col1_header': '',
      'no_beeswarm': False,
      'sortRows': False,
      'format': '{:,.0f}',
      'max_table_rows': 20,
      'save_file': True
    }

    # Add a report section with the table
    self.add_section(
      name = section_name if section_name else '',
      anchor = id + '_anchor',
      description = description if description else '',
      plot = table.plot(data, headers, table_config)
    )
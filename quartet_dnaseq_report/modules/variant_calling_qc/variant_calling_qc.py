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
    ## historical batches' performance
    ref_fn = "history.txt"
    reference_path = os.path.join(os.path.dirname(__file__), 'reference_datasets', ref_fn)
    history_df = pd.read_csv(reference_path, sep='\t')
    if len(history_df) == 0:
      log.debug('No file matched: variant_calling_qc/reference_datasets - history.txt')
    
    ### SUMMARY TABLE 1
    snv_indel_df = pd.DataFrame(); #i=0
    for f in self.find_log_files('variant_calling_qc/snv_indel_summary'):
      if f is None:
        log.debug('No file matched: variant_calling_qc - variants.calling.qc.txt')
      else:
        #i = i+1
        f_p = '%s/%s' % (f['root'], f['fn'])
        tmp_df = pd.read_csv(f_p, sep='\t')
        tmp_df[['SNV precision', 'INDEL precision', 'SNV recall', 'INDEL recall']] = round(tmp_df[['SNV precision', 'INDEL precision', 'SNV recall', 'INDEL recall']]/100, 2)
        snv_indel_df = pd.concat([snv_indel_df, tmp_df], axis=0)
    
    snv_indel_df = snv_indel_df.reset_index(drop=True)
    if len(snv_indel_df) != 0:
      df = snv_indel_df
      # Generate F1-score
      df['SNV F1-score'] = 2*df['SNV precision']*df['SNV recall']/(df['SNV precision']+df['SNV recall'])
      df['INDEL F1-score'] = 2*df['INDEL precision']*df['INDEL recall']/(df['INDEL precision']+df['INDEL recall'])
      # Get the rank of this batch based on F1-score
      snv_rank = self.get_snv_indel_rank(history_df, df, "SNV")
      indel_rank = self.get_snv_indel_rank(history_df, df, "INDEL")
      
      # Generate the mean ± format
      snv_precision = "%s ± %s" % (round(df['SNV precision'].mean(), 2), round(np.std(df['SNV precision']), 2))
      snv_recall = "%s ± %s" % (round(df['SNV recall'].mean(), 2), round(np.std(df['SNV recall']), 2))
      snv_f1 = "%s ± %s" % (round(df['SNV F1-score'].mean(), 2), round(np.std(df['SNV F1-score']), 2))
      indel_precision = "%s ± %s" % (round(df['INDEL precision'].mean(), 2), round(np.std(df['INDEL precision']), 2))
      indel_recall = "%s ± %s" % (round(df['INDEL recall'].mean(), 2), round(np.std(df['INDEL recall']), 2))
      indel_f1 = "%s ± %s" % (round(df['INDEL F1-score'].mean(), 2), round(np.std(df['INDEL F1-score']), 2))

      summary1_df = pd.DataFrame([["SNV", snv_precision, snv_recall, snv_f1, snv_rank], ["INDEL", indel_precision, indel_recall, indel_f1, indel_rank]], columns = ["Type", "Precision", "Recall", "F1-score", "Rank"])
      summary1_path = os.path.join(config.output_dir, "variant_calling_qc_summary.txt")
      summary1_df.to_csv(summary1_path, sep="\t", index=0)
      
      # Transfer into the required format
      summary1_dic = self.convert_input_data_format(summary1_df, 'Type')
      # Plot detailed numbers of performance assessment based on reference datasets
      self.assessment_1('variant_calling_qc_summary', summary1_dic)
      

    ### SUMMARY TABLE 2
    mendelian_df = pd.DataFrame()
    for f in self.find_log_files('variant_calling_qc/mendelian_summary'):
      if f is None:
        log.debug('No file matched: variant_calling_qc - project.summary.txt')
      else:
        f_p = '%s/%s' % (f['root'], f['fn'])
        tmp_df = pd.read_csv(f_p, sep='\t')
        mendelian_df = pd.concat([mendelian_df, tmp_df], axis=0)
    
    mendelian_df = mendelian_df.reset_index(drop=True)
    if len(mendelian_df) != 0:
      df = mendelian_df
      Family = pd.DataFrame(["Family 1", "Family 2", "Family 3"], columns = ["Family"])
      # Extract SNV
      snv_tmp = df[df.Family.str.contains("SNV$")].drop(["Family"], axis = 1).reset_index(drop=True)
      snv_tmp = pd.concat([Family, snv_tmp], axis = 1).dropna(axis=0)
      snv_mcr = "%s ± %s" % (round(snv_tmp['Mendelian_Concordance_Rate'].mean(), 2), round(np.std(snv_tmp['Mendelian_Concordance_Rate']), 2))
      snv_rank = self.get_mendelian_rank(history_df, snv_tmp, "SNV")
      # Extract INDEL
      indel_tmp = df[df.Family.str.contains("INDEL$")].drop(["Family"], axis = 1).reset_index(drop=True)
      indel_tmp = pd.concat([Family, indel_tmp], axis = 1).dropna(axis=0)
      indel_mcr = "%s ± %s" % (round(indel_tmp['Mendelian_Concordance_Rate'].mean(), 2), round(np.std(snv_tmp['Mendelian_Concordance_Rate']), 2))
      indel_rank = self.get_mendelian_rank(history_df, indel_tmp, "INDEL")

      summary2_df = pd.DataFrame([["SNV", snv_mcr, snv_rank], ['INDEL', indel_mcr, indel_rank]], columns=["Type", "Mendelian Concordance Rate", "Rank"])
      summary2_path = os.path.join(config.output_dir, "mendelian_summary.txt")
      summary2_df.to_csv(summary2_path, sep="\t", index=0)

      # Transfer into the required format
      summary2_dic = self.convert_input_data_format(summary2_df, 'Type')
      self.assessment_2('mendelian_summary', summary2_dic)
    

    ### Figure data
    # F1-score of this batch
    f1_snv = snv_indel_df[["Sample", "SNV F1-score"]]
    f1_snv.columns = ["Sample", "F1-score"]
    f1_snv["Type"] = "SNV"
    snv_mcr = []
    for index, row in f1_snv.iterrows():
      if "_1_" in row.Sample:
        MCR = snv_tmp[snv_tmp.Family == "Family 1"]["Mendelian_Concordance_Rate"].to_list()[0]
      elif "_2_" in row.Sample:
        MCR = snv_tmp[snv_tmp.Family == "Family 2"]["Mendelian_Concordance_Rate"].to_list()[0]
      elif "_3_" in row.Sample:
        MCR = snv_tmp[snv_tmp.Family == "Family 3"]["Mendelian_Concordance_Rate"].to_list()[0]
      snv_mcr.append(MCR)
    #print(snv_indel_df); print(f1_snv); print(snv_mcr)
    f1_snv_mcr = pd.concat([f1_snv, pd.DataFrame(snv_mcr, columns = ["MCR"])], axis = 1)
    
    f1_indel = snv_indel_df[["Sample", "INDEL F1-score"]]
    f1_indel.columns = ["Sample", "F1-score"]
    f1_indel["Type"] = "INDEL"
    indel_mcr = []
    for index, row in f1_indel.iterrows():
      if "_1_" in row.Sample:
        MCR = indel_tmp[indel_tmp.Family == "Family 1"]["Mendelian_Concordance_Rate"].to_list()[0]
      elif "_2_" in row.Sample:
        MCR = indel_tmp[indel_tmp.Family == "Family 2"]["Mendelian_Concordance_Rate"].to_list()[0]
      elif "_3_" in row.Sample:
        MCR = indel_tmp[indel_tmp.Family == "Family 3"]["Mendelian_Concordance_Rate"].to_list()[0]
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
    
    ### PLOT 0: png output
    fig_path = os.path.join(config.output_dir, "snv_indel_performance.png")
    sns.set_style("ticks")
    sns.set_context("notebook", font_scale=1.4)
    g = sns.FacetGrid(figure_data, col="Type", hue='Batch')
    g.map_dataframe(sns.scatterplot, x='F1-score', y='MCR', s = 120)
    g.set_titles(col_template="{col_name}")
    g.set_axis_labels("F1-score", "Mendelian Concordance Rate")
    sns.despine(top=False, right=False, left=False, bottom=False)
    handles = g._legend_data.values()
    labels = g._legend_data.keys()
    g.fig.legend(title='', handles=handles, labels=labels, loc='lower center', ncol=2, bbox_to_anchor=(0.47, -0.015, 0, 0),frameon=False) 
    g.fig.subplots_adjust(bottom=0.2)
    g.fig.set_size_inches(9, 5)
    g.savefig(fig_path, dpi = 400)

    ### PLOT 1: SNV Performance
    snv_fig_data = figure_data[figure_data.Type == "SNV"]
    self.quartet_scatter_plot(snv_fig_data, "snv_performance", "SNV", "SNV Performance", "snv_performance")    

    ### PLOT 2: INDEL Performance
    indel_fig_data = figure_data[figure_data.Type == "INDEL"]
    self.quartet_scatter_plot(indel_fig_data, "indel_performance", "INDEL", "INDEL Performance", "INDEL_performance")

    ### DETAILED TABLE 1: snv_indel_summary (variants.calling.qc.txt)
    # Plot detailed numbers of performance assessment based on reference datasets
    # Transfer into the required format    
    detail1_dic = self.convert_input_data_format(snv_indel_df, 'Sample')
    self.detail_1('variant_calling_qc_details', detail1_dic)

    detail1_path = os.path.join(config.output_dir, "variant_calling_qc_details.txt")
    snv_indel_df.to_csv(detail1_path, sep="\t", index=0)

    ### DETAILED TABLE 2: mendelian_summary (mendelian.txt)
    snv_tmp.columns = ["Family", "SNV Detected Variants", "SNV Mendelian Consistent Variants", "SNV MCR Rate"]
    indel_tmp.columns = ["Family", "INDEL Detected Variants", "INDEL Mendelian Consistent Variants", "INDEL MCR Rate"]
    # Combine SNV & INDEL
    mendelian_df_comb = pd.merge(snv_tmp, indel_tmp, on = ["Family"])
    mendelian_df_comb = mendelian_df_comb.drop(['SNV MCR Rate', 'INDEL MCR Rate'], axis=1)
    # Transfer into the required format
    detail2_dic = self.convert_input_data_format(mendelian_df_comb, 'Family')
    self.detail_2('mendelian_details', detail2_dic)
    
    detail2_path = os.path.join(config.output_dir, "mendelian_details.txt")
    mendelian_df_comb.to_csv(detail2_path, sep="\t", index=0)
  
  ### Functions for getting the rank of submitted data
  def get_snv_indel_rank(self, history_df, df, typ):
    history = history_df[history_df.type == typ].f1.to_list()
    history.sort()
    if typ == "SNV":
      query = df['SNV F1-score'].mean()
    else:
      query = df['INDEL F1-score'].mean()
    rank = len(history) + 1
    for i in history:
      if query > i:
        rank = rank - 1
    
    rank_in_history = "%s/%s" % (rank, len(history) + 1)
    return(rank_in_history)

  def get_mendelian_rank(self, history_df, df, typ):
    history = history_df[history_df.type == typ].mendelian.to_list()
    history.sort()
    query = df['Mendelian_Concordance_Rate'].mean()
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
      row = df.iloc[index].tolist()
      convert_df.append(dict(zip(keys, row)))
    
    convert_dic = {}
    for i in convert_df:
      key = i[col_name]
      pop_i = i.pop(col_name)
      convert_dic[key] = i
    
    return(convert_dic)
  

  # Functions for tables and scatter plots
  ## Assessment based on reference datasets (v202103)
  def assessment_1(self, id, data, title='Assessment based on reference datasets', section_name='Assessment based on reference datasets', description="The reference dataset version is v202103.", helptext=None):
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
      'scale': False
    }

    headers['F1-score'] = {
      'title': 'F1-score',
      'description': 'F1-score',
      'scale': False
    }

    headers['Rank'] = {
      'title': 'Rank',
      'description': 'Rank',
      'scale': False
    }

    table_config = {
      'namespace': 'variant_calling_qc_summary',
      'id': id,
      'table_title': 'Assessment based on reference datasets',
      'col1_header': '',
      'no_beeswarm': False,
      'sortRows': False,
      'format': '{:.0f}',
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
  
  ## Assessment based on Quartet family-dependent built-in genetic truth
  def assessment_2(self, id, data, title='Assessment based on Quartet family-dependent built-in genetic truth', section_name='Assessment based on Quartet family-dependent built-in genetic truth', description="", helptext=None):
    """ Create the HTML for assessment based on reference datasets """
    
    headers = OrderedDict()
    headers['Mendelian Concordance Rate'] = {
      'title': 'Mendelian Concordance Rate',
      'description': 'Mendelian Concordance Rate',
      'scale': False
    }

    headers['Rank'] = {
      'title': 'Rank',
      'description': 'Rank',
      'scale': False
    }

    table_config = {
      'namespace': 'mendelian_summary',
      'id': id,
      'table_title': 'Assessment based on Quartet family-dependent built-in genetic truth',
      'col1_header': '',
      'no_beeswarm': False,
      'sortRows': False,
      'format': '{:.0f}',
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
      "use_legend": True
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
  def detail_1(self, id, data, title='Detailed numbers of performance assessment based on reference datasets', section_name='Detailed numbers of performance assessment based on reference datasets', description="", helptext=None):
    """ Create the HTML for detailed numbers of performance assessment based on reference datasets """
    headers = OrderedDict()
    headers['SNV number'] = {
      'title': 'SNV Num',
      'description': 'SNV Total Number',
      'scale': False
    }

    headers['SNV precision'] = {
      'title': 'SNV precision',
      'description': 'SNV Precision',
      'scale': False,
      'format': '{:.2f}'
    }

    headers['SNV recall'] = {
      'title': 'SNV recall',
      'description': 'SNV Recall',
      'scale': False,
      'format': '{:.2f}'
    }

    headers['SNV F1-score'] = {
      'title': 'SNV F1-score',
      'description': 'SNV F1-score',
      'scale': False,
      'format': '{:.2f}'
    }

    headers['INDEL number'] = {
      'title': 'INDEL Num',
      'description': 'INDEL Total Number',
      'scale': False
    }
    
    headers['INDEL precision'] = {
      'title': 'INDEL precision',
      'description': 'INDEL Precision',
      'scale': False,
      'format': '{:.2f}'
    }

    headers['INDEL recall'] = {
      'title': 'INDEL recall',
      'description': 'INDEL Recall',
      'scale': False,
      'format': '{:.2f}'
    }

    headers['INDEL F1-score'] = {
      'title': 'INDEL F1-score',
      'description': 'INDEL F1-score',
      'scale': False,
      'format': '{:.2f}'
    }

    table_config = {
      'namespace': 'variant_calling_qc_details',
      'id': id,
      'table_title': 'Detailed numbers of performance assessment based on reference datasets',
      'col1_header': 'Sample',
      'no_beeswarm': False,
      'sortRows': False,
      'format': '{:.0f}',
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
  

  ## Plot detailed numbers of performance assessment based on Quartet genetic built-in truth
  def detail_2(self, id, data, title='Detailed numbers of performance assessment based on Quartet genetic built-in truth', section_name='Detailed numbers of performance assessment based on Quartet genetic built-in truth', description="", helptext=None):
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
      'table_title': 'Detailed numbers of performance assessment based on Quartet genetic built-in truth',
      'col1_header': 'Family',
      'no_beeswarm': False,
      'sortRows': False,
      'format': '{:.0f}',
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

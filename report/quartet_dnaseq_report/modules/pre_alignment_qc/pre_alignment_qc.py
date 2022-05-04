#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import io
import logging
import os
import zipfile

from multiqc import config
from multiqc.plots import table, linegraph
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
      name='Pre-alignment Quality Control',
      target='Pre-alignment QC',
      #anchor='pre_alignment_qc',
      #href='https://github.com/clinico-omics/quartet-dnaseq-report',
      info=' is an report module to show the data quality before alignment.'
    )

    # Find and load any input files for pre_alignment_qc
    table_summary = []
    for f in self.find_log_files('pre_alignment_qc/summary'):
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
      self.plot_summary_table('pre_alignment_qc_summary', table_summary_dic)
    else:
      log.debug('No file matched: pre_alignment_qc - pre_alignment_qc_summary.txt')

    # Set up class objects to hold parsed data()
    self.general_stats_data = defaultdict(lambda: dict())

    # Find and parse unzipped FastQC reports 
    self.fastqc_data = dict()
    for f in self.find_log_files('pre_alignment_qc/fastqc_data'):
      s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
      self.parse_fastqc_report(f['f'], s_name, f)

    # Find and parse zipped FastQC reports
    for f in self.find_log_files('pre_alignment_qc/fastqc_zip', filecontents=False):
      s_name = f['fn']
      if s_name.endswith('_fastqc.zip'):
        s_name = s_name[:-11]
      # Skip if we already have this report - parsing zip files is slow..
      if s_name in self.fastqc_data.keys():
        log.debug("Skipping '{}' as already parsed '{}'".format(f['fn'], s_name))
        continue
      try:
        fqc_zip = zipfile.ZipFile(os.path.join(f['root'], f['fn']))
      except Exception as e:
        log.warning("Couldn't read '{}' - Bad zip file".format(f['fn']))
        log.debug("Bad zip file error:\n{}".format(e))
        continue
      # FastQC zip files should have just one directory inside, containing report
      d_name = fqc_zip.namelist()[0]
      try:
        with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as fh:
          r_data = fh.read().decode('utf8')
          self.parse_fastqc_report(r_data, s_name, f)
      except KeyError:
        log.warning("Error - can't find fastqc_raw_data.txt in {}".format(f))

    # Filter to strip out ignored sample names
    self.fastqc_data = self.ignore_samples(self.fastqc_data)

    if len(self.fastqc_data) == 0:
      raise UserWarning
    
    # Add to self.css and self.js to be included in template
    self.css = { 'assets/css/multiqc_fastqc.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_fastqc.css') }
    self.js = { 'assets/js/multiqc_fastqc.js' : os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js') }

    # Colours to be used for plotting lines
    self.status_colours = { 'pass': '#5cb85c', 'warn': '#f0ad4e', 'fail': '#d9534f', 'default': '#999' }
    log.info('Found {} reports'.format(len(self.fastqc_data)))
    
    # Now add each section in order
    self.sequence_quality_plot()
    self.gc_content_plot()

  def parse_fastqc_report(self, file_contents, s_name=None, f=None):
    """ Takes contents from a fastq_data.txt file and parses out required
    statistics and data. Returns a dict with keys 'stats' and 'data'.
    Data is for plotting graphs, stats are for top table. """

    if s_name in self.fastqc_data.keys():
      log.debug('Duplicate sample name found! Overwriting: {}'.format(s_name))
    self.add_data_source(f, s_name)
    self.fastqc_data[s_name] = { 'statuses': dict() }

    # Parse the report
    section = None
    s_headers = None
    self.dup_keys = []
    for l in file_contents.splitlines():
      if l == '>>END_MODULE':
        section = None
        s_headers = None
      elif l.startswith('>>'):
        (section, status) = l[2:].split('\t', 1)
        section = section.lower().replace(' ', '_')
        self.fastqc_data[s_name]['statuses'][section] = status
      elif section is not None:
        if l.startswith('#'):
          s_headers = l[1:].split('\t')
          # Special case: Total Deduplicated Percentage header line
          if s_headers[0] == 'Total Deduplicated Percentage':
            self.fastqc_data[s_name]['basic_statistics'].append({
              'measure': 'total_deduplicated_percentage',
              'value': float(s_headers[1])
            })
          else:
            # Special case: Rename dedup header in old versions of FastQC (v10)
            if s_headers[1] == 'Relative count':
              s_headers[1] = 'Percentage of total'
            s_headers = [s.lower().replace(' ', '_') for s in s_headers]
            self.fastqc_data[s_name][section] = list()

        elif s_headers is not None:
          s = l.split('\t')
          row = dict()
          for (i, v) in enumerate(s):
            v.replace('NaN','0')
            try:
              v = float(v)
            except ValueError:
              pass
            row[s_headers[i]] = v
          self.fastqc_data[s_name][section].append(row)
          # Special case - need to remember order of duplication keys
          if section == 'sequence_duplication_levels':
            try:
              self.dup_keys.append(float(s[0]))
            except ValueError:
              self.dup_keys.append(s[0])

    # Tidy up the Basic Stats
    self.fastqc_data[s_name]['basic_statistics'] = {d['measure']: d['value'] for d in self.fastqc_data[s_name]['basic_statistics']}

    # Calculate the average sequence length (Basic Statistics gives a range)
    length_bp = 0
    total_count = 0
    for d in self.fastqc_data[s_name].get('sequence_length_distribution', {}):
      length_bp += d['count'] * self.avg_bp_from_range(d['length'])
      total_count += d['count']
    if total_count > 0:
      self.fastqc_data[s_name]['basic_statistics']['avg_sequence_length'] = length_bp / total_count

  def avg_bp_from_range(self, bp):
    """ Helper function - FastQC often gives base pair ranges (eg. 10-15)
    which are not helpful when plotting. This returns the average from such
    ranges as an int, which is helpful. If not a range, just returns the int """

    try:
      if '-' in bp:
        maxlen = float(bp.split('-',1)[1])
        minlen = float(bp.split('-',1)[0])
        bp = ((maxlen - minlen)/2) + minlen
    except TypeError:
      pass
    return(int(bp))
  
  def get_status_cols(self, section):
    """ Helper function - returns a list of colours according to the FastQC
    status of this module for each sample. """
    colours = dict()
    for s_name in self.fastqc_data:
      status = self.fastqc_data[s_name]['statuses'].get(section, 'default')
      colours[s_name] = self.status_colours[status]
    return colours
  
  def plot_summary_table(self, id, data, title='Summary metrics', section_name='Summary metrics', description=None, helptext=None):
    """ Create the HTML for pre-alignment qc summary """

    headers = OrderedDict()
    headers['%Dup'] = {
      'title': '% Dups',
      'description': '% Duplicate Reads',
      'max': 100,
      'min': 0,
      'scale': "RdYlGn-rev"
    }

    headers['%GC'] = {
      'title': '% GC',
      'description': 'G/C ratio',
      "min": 0,
      'max': 100,
      'scale': "Oranges"
    }

    headers['Total Sequences (million)'] = {
      'title': 'M Seqs',
      'description': 'Total Sequences (million)',
      "min": 0,
      "scale": "YlGnBu",
    }
    
    headers['%Human'] = {
      'title': '% Human',
      'description': '% Human',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%EColi'] = {
      'title': '% EColi',
      'description': '% EColi',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%Adapter'] = {
      'title': '% Adapter',
      'description': '% Adapter',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%Vector'] = {
      'title': '% Vector',
      'description': '% Vector',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%rRNA'] = {
      'title': '% rRNA',
      'description': '% rRNA',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%Virus'] = {
      'title': '% Virus',
      'description': '% Virus',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%Yeast'] = {
      'title': '% Yeast',
      'description': '% Yeast',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%Mitoch'] = {
      'title': '% Mitoch',
      'description': '% Mitoch',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    headers['%No hits'] = {
      'title': '% No hits',
      'description': '% No hits',
      'max': 100,
      'min': 0,
      'scale': 'Reds-rev'
    }

    table_config = {
      'namespace': 'pre_alignment_qc_summary',
      'id': id,
      'table_title': 'Pre-alignment QC Table Summary',
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
  
  def sequence_quality_plot (self):
    """ Create the HTML for the phred quality score plot """

    data = dict()
    for s_name in self.fastqc_data:
      try:
        data[s_name] = {self.avg_bp_from_range(d['base']): d['mean'] for d in self.fastqc_data[s_name]['per_base_sequence_quality']}
      except KeyError:
        pass
    if len(data) == 0:
      log.debug('sequence_quality not found in FastQC reports')
      return None

    pconfig = {
      'id': 'fastqc_per_base_sequence_quality_plot',
      'title': 'FastQC: Mean quality scores',
      'ylab': 'Phred Score',
      'xlab': 'Position (bp)',
      'ymin': 0,
      'xDecimals': False,
      'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
      'colors': self.get_status_cols('per_base_sequence_quality'),
      'yPlotBands': [
        {'from': 28, 'to': 100, 'color': '#c3e6c3'},
        {'from': 20, 'to': 28, 'color': '#e6dcc3'},
        {'from': 0, 'to': 20, 'color': '#e6c3c3'},
      ]
    }
    self.add_section (
      name = 'Sequence quality histograms',
      anchor = 'fastqc_per_base_sequence_quality',
      description = 'The mean quality value across each base position in the read.',
      helptext = '''
      To enable multiple samples to be plotted on the same graph, only the mean quality
      scores are plotted (unlike the box plots seen in FastQC reports).

      Taken from the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html):

      _The y-axis on the graph shows the quality scores. The higher the score, the better
      the base call. The background of the graph divides the y axis into very good quality
      calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
      The quality of calls on most platforms will degrade as the run progresses, so it is
      common to see base calls falling into the orange area towards the end of a read._
      ''',
      plot = linegraph.plot(data, pconfig)
    )

  def gc_content_plot (self):
    """ Create the HTML for the FastQC GC content plot """

    data = dict()
    data_norm = dict()
    for s_name in self.fastqc_data:
      try:
        data[s_name] = {d['gc_content']: d['count'] for d in self.fastqc_data[s_name]['per_sequence_gc_content']}
      except KeyError:
        pass
      else:
        data_norm[s_name] = dict()
        total = sum( [ c for c in data[s_name].values() ] )
        for gc, count in data[s_name].items():
          try:
            data_norm[s_name][gc] = (count / total) * 100
          except ZeroDivisionError:
            data_norm[s_name][gc] = 0
    if len(data) == 0:
      log.debug('per_sequence_gc_content not found in FastQC reports')
      return None

    pconfig = {
      'id': 'fastqc_per_sequence_gc_content_plot',
      'title': 'FastQC: Per sequence GC content',
      'xlab': '% GC',
      'ylab': 'Percentage',
      'ymin': 0,
      'xmax': 100,
      'xmin': 0,
      'format': '{:,.2f}',
      'yDecimals': False,
      'tt_label': '<b>{point.x}% GC</b>: {point.y}',
      'colors': self.get_status_cols('per_sequence_gc_content'),
      'data_labels': [
        {'name': 'Percentages', 'ylab': 'Percentage'},
        {'name': 'Counts', 'ylab': 'Count'}
      ]
    }

    # Try to find and plot a theoretical GC line
    theoretical_gc = None
    theoretical_gc_raw = None
    theoretical_gc_name = None
    for f in self.find_log_files('pre_alignment_qc/fastqc_theoretical_gc'):
      if theoretical_gc_raw is not None:
        log.warning('Multiple FastQC Theoretical GC Content files found, now using {}'.format(f['fn']))
      theoretical_gc_raw = f['f']
      theoretical_gc_name = f['fn']
    if theoretical_gc_raw is None:
      tgc = getattr(config, 'fastqc_config', {}).get('fastqc_theoretical_gc', None)
      if tgc is not None:
        theoretical_gc_name = os.path.basename(tgc)
        tgc_fn = 'fastqc_theoretical_gc_{}.txt'.format(tgc)
        tgc_path = os.path.join(os.path.dirname(__file__), 'fastqc_theoretical_gc', tgc_fn)
        if not os.path.isfile(tgc_path):
          tgc_path = tgc
        try:
          with io.open (tgc_path, 'r', encoding='utf-8') as f:
            theoretical_gc_raw = f.read()
        except IOError:
          log.warning("Couldn't open FastQC Theoretical GC Content file {}".format(tgc_path))
          theoretical_gc_raw = None
    if theoretical_gc_raw is not None:
      theoretical_gc = list()
      for l in theoretical_gc_raw.splitlines():
        if '# FastQC theoretical GC content curve:' in l:
          theoretical_gc_name = l[39:]
        elif not l.startswith('#'):
          s = l.split()
          try:
            theoretical_gc.append([float(s[0]), float(s[1])])
          except (TypeError, IndexError):
            pass

    desc = '''The average GC content of reads. Normal random library typically have a
    roughly normal distribution of GC content.'''
    if theoretical_gc is not None:
      # Calculate the count version of the theoretical data based on the largest data store
      max_total = max([sum (d.values()) for d in data.values() ])
      esconfig = {
        'name': 'Theoretical GC Content',
        'dashStyle': 'Dash',
        'lineWidth': 2,
        'color': '#000000',
        'marker': { 'enabled': False },
        'enableMouseTracking': False,
        'showInLegend': False,
      }
      pconfig['extra_series'] = [ [dict(esconfig)], [dict(esconfig)] ]
      pconfig['extra_series'][0][0]['data'] = theoretical_gc
      pconfig['extra_series'][1][0]['data'] = [ [ d[0], (d[1]/100.0)*max_total ] for d in theoretical_gc ]
      desc = ' **The dashed black line shows theoretical GC content:** `{}`'.format(theoretical_gc_name)

    self.add_section (
      name = 'Per sequence GC content',
      anchor = 'fastqc_per_sequence_gc_content',
      description = desc,
      helptext = '''
      From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html):

      _This module measures the GC content across the whole length of each sequence
      in a file and compares it to a modelled normal distribution of GC content._

      _In a normal random library you would expect to see a roughly normal distribution
      of GC content where the central peak corresponds to the overall GC content of
      the underlying genome. Since we don't know the the GC content of the genome the
      modal GC content is calculated from the observed data and used to build a
      reference distribution._

      _An unusually shaped distribution could indicate a contaminated library or
      some other kinds of biased subset. A normal distribution which is shifted
      indicates some systematic bias which is independent of base position. If there
      is a systematic bias which creates a shifted normal distribution then this won't
      be flagged as an error by the module since it doesn't know what your genome's
      GC content should be._
      ''',
      plot = linegraph.plot([data_norm, data], pconfig)
    )
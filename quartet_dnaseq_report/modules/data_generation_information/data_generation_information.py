#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from quartet_dnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Data Generation Information',
            target='data_generation_information',
            anchor='data_generation_information',
            href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is an report module to show the basic information about the sequencing data.'
        )

        # Find and load any input files for data_generation_information
        for f in self.find_log_files('data_generation_information/information'):
            information = eval(f['f'])
        
        if len(information) != 0:
            self.plot_information('data_generation_information', information)
        else:
            log.debug('No file matched: data_generation_information - data_generation_information.txt')

    def plot_information(self, id, data, title='', section_name='', description=None, helptext=None):
        html_data = ["<dl class='dl-horizontal'>"]
        for k,v in data.items():
            line = "        <dt style='text-align:left;margin-top:1ex'>{}</dt><dd>{}</dd>".format(k,v)
            html_data.append(line)
        html_data.append("    </dl>")

        html = '\n'.join(html_data)
        
        self.add_section(
            name = '',
            anchor = '',
            description = '',
            plot = html
        )
#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
import os
import base64
import logging
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

def read_image(image):
    with open(image, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
        return encoded_string.decode('utf-8')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Supplementary',
            target='supplementary',
            anchor='supplementary',
            href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' is a module to show the additional information about this quality assessment report.'
        )
        
        html = '''
            <!-- Method -->
            <div class='method'>
                <div class='small-12 columns'>
                    <h3 class='section-header black'>Method</h3>
                    <img src="data:image/png;base64,{image}" title='quartet-dna-pipeline' width='70%' height='70%'/>
                    <p>
                        We accepted fastq files, and used sentieon to call germline small variants. [<a class='reference' href='#ref-1'>1</a>] Whole-genome sequencing quality control consists of pre-alignment, post-alignment and variants calling quality control.
                    </p>
                    <p>
                        Pre-alignment quality control focuses on raw fastq files and helps to determine systematic bias and library issue, such as sequencing quality issue, high GC or AT, PCR bias, adapter contaminant, cross species contamination. Fastqc [<a class='reference' href='#ref-2'>2</a>] and fastqscreen [<a class='reference' href='#ref-3'>3</a>] are used to evaluate raw reads quality.
                    </p>
                    <p>
                        Post-alignment quality control focuses on bam files and helps to measure library performance and sample variance, such as sequencing error rate, sequencing depth and coverage consistency. Qualimap [<a class='reference' href='#ref-4'>4</a>] is used to evaluate quality of bam files.
                    </p>
                    <p>
                        Variants calling quality control is to examine accuracy of detected variants based on reference datasets, and estimate potential sequence errors by reproducibility of monozygotic twin daughters and mendelian concordant ratio of Quartet family.
                    </p>
                </div>
            </div>

            <!-- Reference -->
            <div class='reference'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>References</h3>
                <p class='reference'><a name='ref-1'>1.</a> <a href='https://support.sentieon.com/versions/201808.01/manual/'>https://support.sentieon.com/versions/201808.01/manual/</a></p>
                <p class='reference'><a name='ref-2'>2.</a> <a href='https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'>https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></p>
                <p class='reference'><a name='ref-3'>3.</a> <a href='https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/'>https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/</a></p>
                <p class='reference'><a name='ref-4'>4.</a> <a href='http://qualimap.bioinfo.cipf.es/'>http://qualimap.bioinfo.cipf.es/</a></p>
                </div>
            </div>

            <!-- Software version -->
            <div class='software'>
                <h3 class='section-header black'>Software</h3>
                <dl class='dl-horizontal'>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>sentieon (Fastq to VCF)</dt><dd>v2018.08.01</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Fastqc</dt><dd>v0.11.5</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Fastqscreen</dt><dd>v0.12.0</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Qualimap</dt><dd>v2.0.0</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>hap.py (reference datasets benchmark)</dt><dd>v0.3.7</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>VBT (mendelian anaysis)</dt><dd>v1.1</dd>
                </dl>
            </div>

            <!-- Contact us -->
            <div class='contact'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>Contact Us</h3>
                <b>Fudan University Pharmacogenomics Research Center</b>
                <p><strong>Project Manager Ren Luyao</strong></p>
                <li style='margin-top:1ex'>Phone: 15200852771</li>
                <li style='margin-top:1ex'>Email: 18110700050@fudan.edu.cn</li>
                </div>
            </div>

            <!-- Disclaimer -->
            <div class='disclaimer'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>Disclaimer</h3>
                <p>This quality control report is only for this specific test data set and doesn’t represent an evaluation of the business level of the sequencing company. This report is only used for scientific research, not for clinical or commercial use. We don’t bear any economic and legal liabilities for any benefits or losses (direct or indirect) from using the results of this report.</p>
                </div>
            </div>
            '''.format(image=read_image(os.path.join(os.path.dirname(__file__), 'assets', 'img', 'quartet-dna-pipeline_mqc.png')))

        self.add_section(
            name = '',
            anchor = '',
            description = '',
            plot = html
        )
#!/usr/bin/env python

""" Quartet DNAseq Report plugin module """

from __future__ import print_function
import os
import base64
import logging
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

def read_image(image):
    with open(image, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
        return encoded_string.decode('utf-8')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
                
        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None
        
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Supplementary',
            target='The additional information',
            #anchor='supplementary',
            #href='https://github.com/clinico-omics/quartet-dnaseq-report',
            info=' about this quality assessment report.'
        )
        
        html = '''
            <!-- Methods -->
            <div class='methods'>
                <div class='small-12 columns'>
                    <h3 class='section-header black'>Methods</h3>
                    <p>
                        (1)	Tested call sets were compared with benchmark small variants using hap.py (https://github.com/Illumina/hap.py).  Precision is the fraction of called variants in the test dataset that are true, and recall is the fraction of true variants are called in the test dataset. True Positives (TP) are true variants detected in the test dataset. False Negatives (FN) are variants in the reference dataset failed to be detected in the test dataset. False Positive (FP) are variants called in the test dataset but not included in the reference dataset. Precision and recall are defined as below:
                    </p>
                    <p>
                    <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=AM_HTMLorMML-full"></script>
                    <body>
                        <div class="formula">
                            <span>Precision=TP/(TP+FP), Recall=TP/(TP+FN), F1=(2×Precision×Recall)/(Precision+Recall)</span>
                        </div>
                    </body>
                    </p>                  
                    <p>
                        (2)	Mendelian concordance rate (MCR) is the number of variants following Mendelian inheritance laws divided by the total number of variants called among the four Quartet samples. Mendelian concordant variants are the variants shared by the twins (D5 and D6) and following Mendelian inheritance laws with parents (Father: F7 and Mother M8).  Mendelian analysis was performed using VBT (https://github.com/sbg/VBT-TrioAnalysis). When calculating Mendelian concordance rate of small variants, variants on large deletions were not included, because VBT takes these variants as Mendelian violations.
                    </p>
                    <p>
                        (3) Total score = (1+0.5^2) x SNV_score x INDEL_score / (0.5^2 x SNV_score + INDEL_score). SNV_score and INDEL_score are obtained by calculating the mean values of Precision, Recall, and MCR, respectively.
                    </p>
                </div>
            </div>

            <!-- Pipeline -->
            <div class='pipeline'>
                <div class='small-12 columns'>
                    <h3 class='section-header black'>Pipeline</h3>
                    <p>
                        We accepted fastq files, and used Sentieon Genomics to call germline small variants. [<a class='reference' href='#ref-1'>1</a>] The quality control consists of pre-alignment, post-alignment and variants calling quality control. Pre-alignment quality control focuses on raw fastq files and helps to determine systematic bias and library issue, such as sequencing quality issue, high GC or AT, PCR bias, adapter contaminant, cross species contamination. FastQC [<a class='reference' href='#ref-2'>2</a>] and FastQ Screen [<a class='reference' href='#ref-3'>3</a>] are used to evaluate raw reads quality. Post-alignment quality control focuses on bam files and helps to measure library performance and sample variance, such as sequencing error rate, sequencing depth and coverage consistency. Qualimap [<a class='reference' href='#ref-4'>4</a>] is used to evaluate quality of bam files. Variants calling quality control is to examine accuracy of detected variants based on reference datasets, and estimate potential sequence errors by reproducibility of monozygotic twin daughters and mendelian concordant ratio of Quartet family.
                    </p>
                    <img src="data:image/png;base64,{image}" title='quartet-dna-pipeline' width='100%' height='100%'/>
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
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Sentieon Genomics (FASTQ to VCF)</dt><dd>v2019.11.28</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>FastQC</dt><dd>v0.11.5</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>FastQ Screen</dt><dd>v0.12.0</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Qualimap</dt><dd>v2.0.0</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>hap.py (reference datasets benchmark)</dt><dd>v0.3.7</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>VBT (mendelian anaysis)</dt><dd>v1.1</dd>
                </dl>
            </div>

            <!-- Contact us -->
            <div class='contact'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>Contact us</h3>
                <b>Fudan University Pharmacogenomics Research Center</b>
                <li style='margin-top:1ex'>Project manager: Ren Luyao</li>
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
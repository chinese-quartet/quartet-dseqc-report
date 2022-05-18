#!/usr/bin/env python
"""
MultiReport for Quartet DNAseq QC
"""

from setuptools import setup, find_packages

version = '0.2.2'

setup(
    name = 'quartet_dnaseq_report',
    version = version,
    author = 'Yaqing Liu',
    author_email = 'liuyaqing@outlook.com',
    description = 'MultiReport for Quartet DNAseq QC.',
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/chinese-quartet/quartet-dseqc-report',
    download_url = 'https://github.com/chinese-quartet/quartet-dseqc-report/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [
        'multiqc==1.11',
        'plotly==4.9.0',
        'pandas==1.2.4',
        'seaborn==0.11.2',
        'Cython==0.29.28'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'general_information = quartet_dnaseq_report.modules.general_information:MultiqcModule',
            'conclusion = quartet_dnaseq_report.modules.conclusion:MultiqcModule',
            'pre_alignment_qc = quartet_dnaseq_report.modules.pre_alignment_qc:MultiqcModule',
            'post_alignment_qc = quartet_dnaseq_report.modules.post_alignment_qc:MultiqcModule',
            'variant_calling_qc = quartet_dnaseq_report.modules.variant_calling_qc:MultiqcModule',
            'supplementary = quartet_dnaseq_report.modules.supplementary:MultiqcModule'
        ],
        'multiqc.hooks.v1': [
            'execution_start = quartet_dnaseq_report.custom_code:quartet_dnaseq_report_execution_start'
        ],
        'multiqc.cli_options.v1': [
            'disable_plugin = quartet_dnaseq_report.cli:disable_plugin'
        ],
        'multiqc.templates.v1': [
            'report_templates = quartet_dnaseq_report.templates.default'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
#!/usr/bin/env python
""" quartet-dnaseq-report plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks. 
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.quartet_dnaseq_report_version = get_distribution('quartet_dnaseq_report').version


# Add default config options for the things that are used in MultiQC_NGI
def quartet_dnaseq_report_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
        return None
    
    log.info('Running Quartet DNA MultiQC Plugin v{}'.format(config.quartet_dnaseq_report_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    # so we check whether the value is already set. This is to avoid
    # clobbering values that have been customised by users.

    # Module-general_information
    if 'general_information/information' not in config.sp:
        config.update_dict( config.sp, { 'general_information/information': { 'fn_re': r'.*information.json$' } } )
    

    # Module-conclusion
    if 'conclusion/precision_recall_summary' not in config.sp:
        config.update_dict( config.sp, { 'conclusion/precision_recall_summary': { 'fn_re': r'variants.calling.qc.txt$' } } )
    
    if 'conclusion/mendelian_summary' not in config.sp:
        config.update_dict( config.sp, { 'conclusion/mendelian_summary': { 'fn_re': r'.*\.summary.txt$' } } )
    

    # Module-pre_alignment_qc
    if 'pre_alignment_qc/summary' not in config.sp:
        config.update_dict( config.sp, { 'pre_alignment_qc/summary': { 'fn_re': r'pre_alignment.txt$' } } )

    if 'pre_alignment_qc/fastqc_data' not in config.sp:
        config.update_dict( config.sp, { 'pre_alignment_qc/fastqc_data': { 'fn_re': r'fastqc_data.txt$' } } )

    if 'pre_alignment_qc/fastqc_zip' not in config.sp:
        config.update_dict( config.sp, { 'pre_alignment_qc/fastqc_zip': { 'fn_re': r'.*_fastqc.zip' } } )

    if 'pre_alignment_qc/fastqc_theoretical_gc' not in config.sp:
        config.update_dict( config.sp, { 'pre_alignment_qc/fastqc_theoretical_gc': { 'fn_re': r'fastqc_theoretical_gc_hg38_genome.txt$' } } )
    

    # Module-post_alignment_qc
    if 'post_alignment_qc/summary' not in config.sp:
        config.update_dict( config.sp, { 'post_alignment_qc/summary': { 'fn_re': r'post_alignment.txt$' } } )

    if 'post_alignment_qc/bamqc/genome_results' not in config.sp:
        config.update_dict( config.sp, { 'post_alignment_qc/bamqc/genome_results': { 'fn_re': r'^genome_results.txt$' } } )

    if 'post_alignment_qc/bamqc/coverage' not in config.sp:
        config.update_dict( config.sp, { 'post_alignment_qc/bamqc/coverage': { 'fn_re': r'coverage_histogram.txt$' } } )

    if 'post_alignment_qc/bamqc/insert_size' not in config.sp:
        config.update_dict( config.sp, { 'post_alignment_qc/bamqc/insert_size': { 'fn_re': r'insert_size_histogram.txt$' } } )

    if 'post_alignment_qc/bamqc/genome_fraction' not in config.sp:
        config.update_dict( config.sp, { 'post_alignment_qc/bamqc/genome_fraction': { 'fn_re': r'genome_fraction_coverage.txt$' } } )

    if 'post_alignment_qc/bamqc/gc_dist' not in config.sp:
        config.update_dict( config.sp, { 'post_alignment_qc/bamqc/gc_dist': { 'fn_re': r'mapped_reads_gc-content_distribution.txt$' } } )
    

    # Module-variant_calling_qc
    if 'variant_calling_qc/precision_recall_summary' not in config.sp:
        config.update_dict( config.sp, { 'variant_calling_qc/precision_recall_summary': { 'fn_re': r'variants.calling.qc.txt$' } } )
    
    if 'variant_calling_qc/mendelian_summary' not in config.sp:
        config.update_dict( config.sp, { 'variant_calling_qc/mendelian_summary': { 'fn_re': r'.*\.summary.txt$' } } )
    
    
    config.module_order = ['general_information', 'conclusion', 'pre_alignment_qc', 'post_alignment_qc', 'variant_calling_qc', 'supplementary']

    config.exclude_modules = ['fastqc', 'fastq_screen', 'qualimap']
    
    config.log_filesize_limit = 2000000000
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

    log.info('Running Example MultiQC Plugin v{}'.format(config.quartet_dnaseq_report_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    # so we check whether the value is already set. This is to avoid
    # clobbering values that have been customised by users.

    # Add to the search patterns used by data_generation_information
    if 'data_generation_information/information' not in config.sp:
        config.update_dict( config.sp, { 'data_generation_information/information': { 'fn_re': '^data_generation_information.json$' } } )


    # Add to the search patterns used by performance_assessment
    if 'performance_assessment/precision_recall_snv' not in config.sp:
        config.update_dict( config.sp, { 'performance_assessment/precision_recall_snv': { 'fn_re': '^precision_recall_snv.txt$' } } )

    if 'performance_assessment/precision_recall_indel' not in config.sp:
        config.update_dict( config.sp, { 'performance_assessment/precision_recall_indel': { 'fn_re': '^precision_recall_indel.txt$' } } )

    if 'performance_assessment/mendelian_jaccard_index_snv' not in config.sp:
        config.update_dict( config.sp, { 'performance_assessment/mendelian_jaccard_index_snv': { 'fn_re': '^mendelian_jaccard_index_snv.txt$' } } )

    if 'performance_assessment/mendelian_jaccard_index_indel' not in config.sp:
        config.update_dict( config.sp, { 'performance_assessment/mendelian_jaccard_index_indel': { 'fn_re': '^mendelian_jaccard_index_indel.txt$' } } )
    
    
    # Add to the search patterns used by prealignment_qc
    if 'prealignment_qc/prealignment_qc_summary' not in config.sp:
        config.update_dict( config.sp, { 'prealignment_qc/prealignment_qc_summary': { 'fn_re': '^prealignment_qc_summary.txt$' } } )

    if 'prealignment_qc/fastqc_data' not in config.sp:
        config.update_dict( config.sp, { 'prealignment_qc/fastqc_data': { 'fn_re': 'fastqc_data.txt' } } )

    if 'prealignment_qc/fastqc_zip' not in config.sp:
        config.update_dict( config.sp, { 'prealignment_qc/fastqc_zip': { 'fn_re': '.*_fastqc.zip' } } )

    if 'prealignment_qc/fastqc_theoretical_gc' not in config.sp:
        config.update_dict( config.sp, { 'prealignment_qc/fastqc_theoretical_gc': { 'fn_re': '.*fastqc_theoretical_gc.*' } } )

    # Add to the search patterns used by postalignment_qc
    if 'postalignment_qc/postlignment_qc_summary' not in config.sp:
        config.update_dict( config.sp, { 'postalignment_qc/postalignment_qc_summary': { 'fn_re': '^postalignment_qc_summary.txt$' } } )

    if 'postalignment_qc/bamqc/genome_results' not in config.sp:
        config.update_dict( config.sp, { 'postalignment_qc/bamqc/genome_results': { 'fn_re': '^genome_results.txt$' } } )

    if 'postalignment_qc/bamqc/coverage' not in config.sp:
        config.update_dict( config.sp, { 'postalignment_qc/bamqc/coverage': { 'fn_re': '^coverage_histogram.txt$' } } )

    if 'postalignment_qc/bamqc/insert_size' not in config.sp:
        config.update_dict( config.sp, { 'postalignment_qc/bamqc/insert_size': { 'fn_re': '^insert_size_histogram.txt$' } } )

    if 'postalignment_qc/bamqc/genome_fraction' not in config.sp:
        config.update_dict( config.sp, { 'postalignment_qc/bamqc/genome_fraction': { 'fn_re': '^genome_fraction_coverage.txt$' } } )

    if 'postalignment_qc/bamqc/gc_dist' not in config.sp:
        config.update_dict( config.sp, { 'postalignment_qc/bamqc/gc_dist': { 'fn_re': '^mapped_reads_gc-content_distribution.txt$' } } )

    # Add to the search patterns used by variant_calling_qc
    if 'variant_calling_qc/variant_calling_qc_summary' not in config.sp:
        config.update_dict( config.sp, { 'variant_calling_qc/variant_calling_qc_summary': { 'fn_re': '^variant_calling_qc_summary.txt$' } } )

    if 'variant_calling_qc/benchmark_score' not in config.sp:
        config.update_dict( config.sp, { 'variant_calling_qc/benchmark_score': { 'fn_re': '^benchmark_score.txt$' } } )
    
    # # Some additional filename cleaning
    # config.fn_clean_exts.extend([
    #     '.my_tool_extension',
    #     '.removeMetoo'
    # ])

    # # Ignore some files generated by the custom pipeline
    # config.fn_ignore_paths.extend([
    #     '*/my_awesome_pipeline/fake_news/*',
    #     '*/my_awesome_pipeline/red_herrings/*',
    #     '*/my_awesome_pipeline/noisy_data/*',
    #     '*/my_awesome_pipeline/rubbish/*'
    # ])

    config.module_order = ['data_generation_information', 'performance_assessment', 'prealignment_qc', 'postalignment_qc', 'variant_calling_qc', 'supplementary']

    config.log_filesize_limit = 2000000000
#!/usr/bin/env python3

import os
import re
import json
import click
# You may need to install https://github.com/yjcyxky/biominer-app-util firstly.
from biominer_app_util.cli import render_app
from subprocess import Popen, PIPE


def read_json(json_file):
    with open(json_file, "r") as f:
        return json.load(f)


def write_json(data, json_file):
    with open(json_file, 'w') as f:
        json.dump(data, f)


@click.group()
def dseqc():
    pass


@dseqc.command(help="Run the pipeline for DNA-Seq (WGS/WES) data. You need to specify the --bed-file argument if you want to analyze WES data.")
@click.option('--d5-r1', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="D5 Read1 File.")
@click.option('--d5-r2', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="D5 Read2 File.")
@click.option('--d6-r1', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="D6 Read1 File.")
@click.option('--d6-r2', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="D6 Read2 File.")
@click.option('--f7-r1', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="F7 Read1 File.")
@click.option('--f7-r2', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="F7 Read2 File.")
@click.option('--m8-r1', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="M8 Read1 File.")
@click.option('--m8-r2', required=True, multiple=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="M8 Read2 File.")
@click.option('--platform', '-p', required=False,
              type=click.Choice(["BGI", "ILLUMINA"]),
              help="Which platform?.", default="ILLUMINA")
@click.option('--bed-file', '-b', required=False,
              type=click.Path(exists=True, dir_okay=False, file_okay=True),
              help="A bed file for your wes data.")
@click.option('--benchmarking-dir', '-B', required=True,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="A directory which contains reference datasets for benchmarking.")
@click.option('--reference-data-dir', '-R', required=True,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="A directory which contains reference data files.")
@click.option('--fastq-screen-dir', '-F', required=True,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="A directory which contains fastq_screen reference files.")
@click.option('--sentieon-server', '-S', required=True,
              help="A url for sentieon license server.")
@click.option('--output-dir', required=False,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="The output directory.")
def fq_workflow(d5_r1, d5_r2, d6_r1, d6_r2, f7_r1, f7_r2, m8_r1, m8_r2,
                platform, bed_file, output_dir, benchmarking_dir, fastq_screen_dir, reference_data_dir, sentieon_server):
    for item in [d5_r1, d6_r1, f7_r1, m8_r1]:
        if not re.match(r'.*_R1.(fastq|fq).gz', item):
            raise Exception(
                "The file (%s) must be with suffixes of _R1.fastq.gz or _R1.fq.gz" % item)

    for item in [d5_r2, d6_r2, f7_r2, m8_r2]:
        if not re.match(r'.*_R2.(fastq|fq).gz', item):
            raise Exception(
                "The file (%s) must be with suffixes of _R2.fastq.gz or _R2.fq.gz" % item)

    if bed_file:
        wdl_dir = '/venv/wes-workflow'
    else:
        wdl_dir = '/venv/wgs-workflow'

    if not os.path.exists(wdl_dir):
        print("Cannot find the workflow, please contact the administrator.")

    project_name = "dseqc"
    data_dict = {
        "project_name": project_name,
        "pl": platform,
        "fastq_or_vcf": "fastq",
        "benchmarking_dir": benchmarking_dir,
        "benchmarking_region": os.path.join(benchmarking_dir, "Quartet.high.confidence.region.v202103.bed"),
        "screen_ref_dir": fastq_screen_dir,
        "dbsnp_dir": reference_data_dir,
        "reference_bed_dict": os.path.join(reference_data_dir, "GRCh38.d1.vd1.dict"),
        "dbmills_dir": reference_data_dir,
        "ref_dir": reference_data_dir,
        "SENTIEON_LICENSE": sentieon_server,
        "fastq_1_D5": d5_r1,
        "fastq_2_D5": d5_r2,
        "fastq_1_D6": d6_r1,
        "fastq_2_D6": d6_r2,
        "fastq_1_F7": f7_r1,
        "fastq_2_F7": f7_r2,
        "fastq_1_M8": m8_r1,
        "fastq_2_M8": m8_r2
    }

    if bed_file:
        data_dict["bed"] = bed_file

    output_workflow_dir = os.path.join(output_dir, project_name)
    os.makedirs(output_workflow_dir, exist_ok=True)

    render_app(wdl_dir, output_dir=output_workflow_dir,
               project_name=project_name, sample=data_dict)

    def call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path):
        cmd = ['java', '-Dconfig.file=/venv/cromwell-local.conf', '-jar',
               '/venv/share/cromwell/cromwell.jar', 'run', workflow_fpath, "-i", inputs_fpath,
               "-p", tasks_path, "--workflow-root", workflow_root]
        print('Run workflow and output results to %s.' % workflow_root)
        proc = Popen(cmd, stdin=PIPE)
        proc.communicate()

    inputs_fpath = os.path.join(output_workflow_dir, "inputs")
    workflow_fpath = os.path.join(output_workflow_dir, "workflow.wdl")
    workflow_root = output_dir
    tasks_path = os.path.join(output_workflow_dir, "tasks.zip")
    call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path)


@dseqc.command(help="Run the pipeline for DNA-Seq data.")
@click.option('--vcf-d5', required=True, multiple=False,
              type=click.Path(exists=True, file_okay=True),
              help="D5 VCF Files.")
@click.option('--vcf-d6', required=True, multiple=False,
              type=click.Path(exists=True, file_okay=True),
              help="D6 VCF Files.")
@click.option('--vcf-f7', required=True, multiple=False,
              type=click.Path(exists=True, file_okay=True),
              help="F7 VCF Files.")
@click.option('--vcf-m8', required=True, multiple=False,
              type=click.Path(exists=True, file_okay=True),
              help="M8 VCF Files.")
@click.option('--platform', '-p', required=False,
              type=click.Choice(["BGI", "ILLUMINA"]),
              help="Which platform?.", default="ILLUMINA")
@click.option('--bed-file', '-b', required=False,
              type=click.Path(exists=True, file_okay=True),
              help="A bed file for your wes data.")
@click.option('--benchmarking-dir', '-B', required=True,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="A directory which contains reference datasets for benchmarking.")
@click.option('--reference-data-dir', '-R', required=True,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="A directory which contains reference data files.")
@click.option('--output-dir', required=False,
              type=click.Path(exists=True, dir_okay=True),
              help="The output directory.")
def vcf_workflow(vcf_d5, vcf_d6, vcf_f7, vcf_m8, platform, bed_file, output_dir, benchmarking_dir, reference_data_dir):
    for item in [vcf_d5, vcf_d6, vcf_f7, vcf_m8]:
        if not re.match(r'.*.vcf', item):
            raise Exception(
                "The file (%s) must be with suffixes of .vcf" % item)

    if bed_file:
        wdl_dir = '/venv/wes-workflow'
    else:
        wdl_dir = '/venv/wgs-workflow'

    if not os.path.exists(wdl_dir):
        print("Cannot find the workflow, please contact the administrator.")

    project_name = "dseqc"
    data_dict = {
        "project_name": project_name,
        "platform": platform,
        "fastq_or_vcf": "vcf",
        "benchmarking_dir": benchmarking_dir,
        "benchmarking_region": os.path.join(benchmarking_dir, "Quartet.high.confidence.region.v202103.bed"),
        "dbsnp_dir": reference_data_dir,
        "reference_bed_dict": os.path.join(reference_data_dir, "GRCh38.d1.vd1.dict"),
        "dbmills_dir": reference_data_dir,
        "ref_dir": reference_data_dir,
        "vcf_D5": vcf_d5,
        "vcf_D6": vcf_d6,
        "vcf_F7": vcf_f7,
        "vcf_M8": vcf_m8,
    }

    if bed_file:
        data_dict["bed"] = bed_file

    output_workflow_dir = os.path.join(output_dir, project_name)
    os.makedirs(output_workflow_dir, exist_ok=True)

    render_app(wdl_dir, output_dir=output_workflow_dir,
               project_name=project_name, sample=data_dict)

    def call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path):
        cmd = ['java', '-Dconfig.file=/venv/cromwell-local.conf',
               '-jar', '/venv/share/cromwell/cromwell.jar', 'run', workflow_fpath, "-i", inputs_fpath,
               "-p", tasks_path, "--workflow-root", workflow_root]
        print('Run workflow and output results to %s.' % workflow_root)
        proc = Popen(cmd, stdin=PIPE)
        proc.communicate()

    inputs_fpath = os.path.join(output_workflow_dir, "inputs")
    workflow_fpath = os.path.join(output_workflow_dir, "workflow.wdl")
    workflow_root = output_dir
    tasks_path = os.path.join(output_workflow_dir, "tasks.zip")
    call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path)


@dseqc.command(help="Run the report for DNA-Seq results.")
@click.option('--result-dir', '-d', required=True,
              type=click.Path(exists=True, file_okay=True),
              help="A directory which contains a series of results from DNA-Seq pipeline.")
@click.option('--output-dir', '-o', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="A directory which will store the output report.")
def report(result_dir, output_dir):
    cmd = ['quartet-dseqc-report', '-d', result_dir, "-o", output_dir]
    print('Run quartet-dseqc-report and output the report to %s.' % output_dir)
    proc = Popen(cmd, stdin=PIPE)
    proc.communicate()


if __name__ == '__main__':
    dseqc()

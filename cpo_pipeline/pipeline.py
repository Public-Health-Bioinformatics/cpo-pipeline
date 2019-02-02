#!/usr/bin/env python

import os
import argparse
import configparser
import csv
import drmaa
import subprocess
import sys
import multiprocessing
from pkg_resources import resource_filename
from pprint import pprint
import cpo_pipeline

def prepare_job(job, session):
    job_template = session.createJobTemplate()
    job_template.jobName = job['job_name']
    job_template.nativeSpecification = job['native_specification']
    job_template.jobEnvironment = os.environ
    job_template.workingDirectory = os.getcwd()
    job_template.remoteCommand = job['remote_command']
    job_template.args = job['args']
    job_template.joinFiles = True
    
    return job_template

def run_jobs(jobs):
    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

def collect_final_outputs(outdir, sample_id):
    final_outputs = {}
    final_outputs['sample_id'] = sample_id
    total_bp = cpo_pipeline.assembly.parsers.result_parsers.parse_total_bp(
        "/".join([
            outdir,
            sample_id,
            'pre-assembly_qc',
            'totalbp'
        ])
    )
    final_outputs['bp'] = total_bp

    [mlst_result] = cpo_pipeline.typing.parsers.result_parsers.parse_mlst_result(
        "/".join([
            outdir,
            sample_id,
            'typing',
            'mlst',
            'mlst.tsv'
        ])
    )
    final_outputs['MLST_SCHEME'] = mlst_result['scheme_id']
    final_outputs['MLST'] = mlst_result['sequence_type']
    allele_number = 1
    for key, value in mlst_result['multi_locus_alleles'].items():
        final_outputs['MLST_ALLELE_' + str(allele_number)] = key + "(" + value + ")"
        allele_number += 1

    return [final_outputs]

def main(args):
    """
    """
    
    config = configparser.ConfigParser()
    config.read(args.config_file)

    plasmids_command_line = [
        'cpo-pipeline',
        'plasmids',
        '--ID', args.sample_id,
        '--R1', args.reads1_fastq,
        '--R2', args.reads2_fastq,
        '--outdir', args.outdir,
    ]

    subprocess.run(assembly_command_line)

    assembly_command_line = [
        'cpo-pipeline',
        'assembly',
        '--ID', args.sample_id,
        '--R1', args.reads1_fastq,
        '--R2', args.reads2_fastq,
        '--outdir', args.outdir,
    ]
    
    subprocess.run(assembly_command_line)
        
    typing_command_line = [
        'cpo-pipeline',
        'typing',
        '--ID', args.sample_id,
        '--assembly', "/".join([args.outdir, args.sample_id, 'assembly', 'contigs.fa']),
        '--outdir', args.outdir,
    ]

    subprocess.run(typing_command_line)

    resistance_command_line = [
        'cpo-pipeline',
        'resistance',
        '--ID', args.sample_id,
        '--assembly', "/".join([args.outdir, args.sample_id, 'assembly', 'contigs.fa']),
        '--outdir', args.outdir,
    ]

    subprocess.run(resistance_command_line)

    final_outputs = collect_final_outputs(args.outdir, args.sample_id)
    pprint(final_outputs)
    final_output_csv_path = "/".join([
        args.outdir,
        args.sample_id,
        'final_output.csv'
    ])

    final_outputs_headers = [
        'sample_id',
        'bp',
        'MLST_SCHEME',
        'MLST',
        'MLST_ALLELE_1',
        'MLST_ALLELE_2',
        'MLST_ALLELE_3',
        'MLST_ALLELE_4',
        'MLST_ALLELE_5',
        'MLST_ALLELE_6',
        'MLST_ALLELE_7',
    ]

    with open(final_output_csv_path, 'w+') as f:
        writer = csv.DictWriter(f, fieldnames=final_outputs_headers, delimiter='\t')
        writer.writeheader()
        for row in final_outputs:
            writer.writerow(row)


def multi(args):
    """
    """
    
    config = configparser.ConfigParser()
    config.read(args.config_file)

    pool = multiprocessing.Pool(args.parallel)

    script_name = os.path.basename(os.path.realpath(sys.argv[0]))

    arguments = []
    with open(args.input_file) as input_file:
        fieldnames = ['sample_id', 'reads1_fastq', 'reads2_fastq']
        reader = csv.DictReader(
            (row for row in input_file if not row.startswith('#')),
            delimiter='\t',
            fieldnames=fieldnames
        )
        for row in reader:
            parser = argparse.ArgumentParser(prog=script_name, description='')
            parser.add_argument("-i", "--ID", dest="sample_id",
                                help="identifier of the isolate")
            parser.add_argument("-1", "--R1", dest="reads1_fastq",
                                help="absolute file path forward read (R1)")
            parser.add_argument("-2", "--R2", dest="reads2_fastq",
                                help="absolute file path to reverse read (R2)")
            parser.add_argument("-o", "--outdir", dest="outdir", default='./',
                                help="absolute path to output folder")
            parser.add_argument('-c', '--config', dest='config_file',
                                default=resource_filename('data', 'config.ini'),
                                help='Config File', required=False)
            job_args = parser.parse_args([
                '--ID', row['sample_id'],
                '--R1', row['reads1_fastq'],
                '--R2', row['reads2_fastq'],
                '--outdir', args.outdir,
            ])
            arguments.append(job_args)

    pool.map(cpo_pipeline.pipeline.main, arguments)
    
if __name__ == '__main__':
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name, description='')
    parser.add_argument("-i", "--ID", dest="sample_id",
                        help="identifier of the isolate", required=True)
    parser.add_argument("-1", "--R1", dest="reads1_fastq",
                        help="absolute file path forward read (R1)", required=True)
    parser.add_argument("-2", "--R2", dest="reads2_fastq",
                        help="absolute file path to reverse read (R2)", required=True)
    parser.add_argument("-o", "--outdir", dest="outdir", default='./',
                        help="absolute path to output folder")
    parser.add_argument('-c', '--config', dest='config_file',
                        default=resource_filename('data', 'config.ini'),
                        help='Config File', required=False)
    args = parser.parse_args()
    main(args)

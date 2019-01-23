#!/usr/bin/env python

import os
import argparse
import configparser
import csv
import subprocess
import sys
import multiprocessing
from pkg_resources import resource_filename

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


def main(args):
    """
    """
    
    config = configparser.ConfigParser()
    config.read(args.config_file)

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
        '--assembly', "/".join(args.outdir, args.sample_id, 'assembly', 'contigs.fa'),
        '--outdir', args.outdir,
    ]

    subprocess.run(typing_command_line)

    resistance_command_line = [
        'cpo-pipeline',
        'resistance',
        '--ID', args.sample_id,
        '--assembly', "/".join(args.outdir, args.sample_id, 'assembly', 'contigs.fa'),
        '--outdir', args.outdir,
    ]

    subprocess.run(resistance_command_line)

def multi(args):
    """
    """
    
    config = configparser.ConfigParser()
    config.read(args.config_file)

    pool = multiprocessing.Pool(8)

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

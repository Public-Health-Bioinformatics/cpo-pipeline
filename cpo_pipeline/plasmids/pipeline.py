#!/usr/bin/env python

'''

Example usage:

pipeline.py -i BC11-Kpn005_S2 --R1 BC11-Kpn005_S2_L001_R1_001.fastq.gz \
  --R2 BC11-Kpn005_S2_L001_R2_001.fastq.gz -o output_dir
'''

import argparse
import configparser
import datetime
import errno
import glob
import os
import re
import sys
import urllib.request

import drmaa

from pkg_resources import resource_filename

from cpo_pipeline.pipeline import prepare_job


def main(args):
    """
    main entrypoint
    Args:
        args():
    Returns:
        (void)
    """

    config = configparser.ConfigParser()
    config.read(args.config_file)

    if args.mash_refseq_plasmid_db and not config['databases']['mash_refseq_plasmid_db']:
        mash_refseq_plasmid_db = args.mash_refseq_plasmid_db
    else:
        mash_refseq_plasmid_db = config['databases']['mash_refseq_plasmid_db']
    if args.mash_custom_plasmid_db and not config['databases']['mash_custom_plasmid_db']:
        mash_custom_plasmid_db = args.mash_custom_plasmid_db
    else:
        mash_custom_plasmid_db = config['databases']['mash_custom_plasmid_db']

    sample_id = args.sample_id
    reads1_fastq = args.reads1_fastq
    reads2_fastq = args.reads2_fastq
    output_dir = args.outdir

    job_script_path = resource_filename('data', 'job_scripts')

    file_paths = {
        "mash_refseq_plasmid_path": "/".join([output_dir, sample_id, "plasmids", "mashscreen.refseq.plasmid.tsv"]),
        "mash_custom_plasmid_path": "/".join([output_dir, sample_id, "plasmids", "mashscreen.custom.plasmid.tsv"]),
    }
    
    mash_jobs = [
        {
            'job_name': 'mash_screen_refseq_plasmid',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_screen.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--queries", mash_refseq_plasmid_db,
                "--output_file", file_paths['mash_refseq_plasmid_path']
            ],
        },
        {
            'job_name': 'mash_screen_custom_plasmid',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_screen_custom_db.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--plasmid-db-dir", mash_custom_plasmid_db,
                "--output_file", file_paths['mash_custom_plasmid_path']
            ],
        },
    ]

    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in mash_jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    # TODO:
    # determine which mash hits are 'candidates'
    # copy the candidate plasmid match fasta files over to working directory
    # index plasmid reference with samtools & bwa
    # align reads to referece with bwa
    # determine depth of alignment with samtools
    # call snps with freebayes
    
    bwa_jobs = [
        {
            'job_name': 'samtools_faidx',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'samtools_faidx.sh'),
            'args': [
                "--fasta", plasmid_fasta,
            ],
        },
    ]
        
if __name__ == "__main__":
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name, description='')
    parser.add_argument("-i", "--ID", dest="sample_id",
                        help="identifier of the isolate", required=True)
    parser.add_argument("-1", "--R1", dest="reads1_fastq",
                        help="absolute file path forward read (R1)", required=True)
    parser.add_argument("-2", "--R2", dest="reads2_fastq",
                        help="absolute file path to reverse read (R2)", required=True)
    parser.add_argument("-o", "--outdir", dest="output", default='./',
                        help="absolute path to output folder")
    parser.add_argument("--mash-refseq-plasmiddb", dest="mash_refseq_plasmid_db",
                        help="absolute path to mash reference database")
    parser.add_argument("--mash-custom-plasmiddb", dest="mash_custom_plasmid_db",
                        help="absolute file path to directory of custom plasmid mash sketches")
    parser.add_argument('-c', '--config', dest='config_file',
                        default=resource_filename('data', 'config.ini'),
                        help='Config File', required=False)
    args = parser.parse_args()
    main(args)

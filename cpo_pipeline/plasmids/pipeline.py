#!/usr/bin/env python

'''

Example usage:

pipeline.py -i BC11-Kpn005_S2 --R1 BC11-Kpn005_S2_L001_R1_001.fastq.gz \
  --R2 BC11-Kpn005_S2_L001_R2_001.fastq.gz -o output_dir
'''

import argparse
import configparser
import csv
import datetime
import errno
import glob
import os
import operator
import re
import shutil
import sys
import urllib.request
from pprint import pprint
import drmaa

from pkg_resources import resource_filename

from cpo_pipeline.pipeline import prepare_job, run_jobs
from cpo_pipeline.assembly.parsers import result_parsers
from cpo_pipeline.plasmids import parsers


def samtools_discrete_jobs(candidates, file_paths):
    samtools_view_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".sam",
        ])
        samtools_view_job = {
            'job_name': 'samtools_view',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_view.sh'),
            # '--f 1540' excludes the following reads:
            # - read unmapped (0x4)
            # - read fails platform/vendor quality checks (0x200)
            # - read is PCR or optical duplicate (0x400)
            'args': [
                "--input", alignment,
                "--flags", 1540,
                "--output", re.sub("\.sam$", ".mapped.dedup.bam", alignment),
            ]
        }
        samtools_view_jobs.append(samtools_view_job)

    run_jobs(samtools_view_jobs)

    samtools_sort_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".mapped.dedup.bam",
        ])
        samtools_sort_job = {
            'job_name': 'samtools_sort',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_sort.sh'),
            'args': [
                "--input", alignment,
                "--name-order",
                "--output", re.sub("\.bam$", ".namesort.bam", alignment),
            ]
        }
        samtools_sort_jobs.append(samtools_sort_job)

    run_jobs(samtools_sort_jobs)

    samtools_fixmate_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".mapped.dedup.namesort.bam",
        ])
        samtools_fixmate_job = {
            'job_name': 'samtools_fixmate',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_fixmate.sh'),
            'args': [
                "--input", alignment,
                "--output", re.sub("\.bam$", ".fixmate.bam", alignment),
            ]
        }
        samtools_fixmate_jobs.append(samtools_fixmate_job)

    run_jobs(samtools_fixmate_jobs)

    samtools_sort_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".mapped.dedup.namesort.fixmate.bam",
        ])
        samtools_sort_job = {
            'job_name': 'samtools_sort',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_sort.sh'),
            'args': [
                "--input", alignment,
                "--output", re.sub("\.bam$", ".coordsort.bam", alignment),
            ]
        }
        samtools_sort_jobs.append(samtools_sort_job)

    run_jobs(samtools_sort_jobs)

    samtools_markdup_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".mapped.dedup.namesort.fixmate.coordsort.bam",
        ])
        samtools_markdup_job = {
            'job_name': 'samtools_markdup',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_markdup.sh'),
            'args': [
                "--input", alignment,
                "--output", re.sub("\.bam$", ".markdup.bam", alignment),
            ]
        }
        samtools_markdup_jobs.append(samtools_markdup_job)

    run_jobs(samtools_markdup_jobs)

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
        "plasmid_output_path": "/".join([output_dir, sample_id, "plasmids"]),
        "mash_refseq_plasmid_path": "/".join([output_dir, sample_id, "plasmids", "mashscreen.refseq.plasmid.tsv"]),
        "mash_custom_plasmid_path": "/".join([output_dir, sample_id, "plasmids", "mashscreen.custom.plasmid.tsv"]),
        "mash_custom_plasmid_candidates_path": "/".join([output_dir, sample_id, "plasmids", "candidates.tsv"]),
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
            'native_specification': '-pe smp 8 -shell y',
            'remote_command': os.path.join(job_script_path, 'mash_screen_custom_db.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--plasmid-db-dir", "/".join([mash_custom_plasmid_db, "mash"]),
                "--output_file", file_paths['mash_custom_plasmid_path']
            ],
        },
    ]

    run_jobs(mash_jobs)

    # TODO:
    # DONE determine which mash hits are 'candidates'
    # DONE copy the candidate plasmid match fasta files over to working directory
    # DONE index plasmid reference with samtools & bwa
    # DONE align reads to referece with bwa
    # DONE determine depth of alignment with samtools
    # call snps with freebayes

    mash_screen_results = result_parsers.parse_mash_result(file_paths['mash_custom_plasmid_path'])
    # pprint(mash_screen_results)

    custom_plasmid_db_data = {}
    for dat_file in glob.glob("/".join([mash_custom_plasmid_db, "data", "*.dat"])):
        [dat] = parsers.custom_plasmid_db_dat_parser(dat_file)
        custom_plasmid_db_data[dat['accession']] = dat

    # pprint(custom_plasmid_db_data)
    for mash_screen_result in mash_screen_results:
        accession = re.sub('\.fna$', '', mash_screen_result['query_id'])
        mash_screen_result['accession'] = accession
        mash_screen_result['allele'] = custom_plasmid_db_data[accession]['allele']
        mash_screen_result['circularity'] = custom_plasmid_db_data[accession]['circularity']
        mash_screen_result['plasmid_length'] = custom_plasmid_db_data[accession]['plasmid_length']
        mash_screen_result['incompatibility_group'] = custom_plasmid_db_data[accession]['incompatibility_group']

    mash_screen_results.sort(key=operator.itemgetter('accession'))
    mash_screen_results.sort(key=operator.itemgetter('plasmid_length'), reverse=True)
    mash_screen_results.sort(key=operator.itemgetter('identity'), reverse=True)
    mash_screen_results.sort(key=operator.itemgetter('circularity'))
    mash_screen_results.sort(key=operator.itemgetter('incompatibility_group'))

    # print("mash_identity\taccession\tcircularity\tplasmid_length\tincompatibility_group")
    candidates_keys = [
        'identity',
        'accession',
        'circularity',
        'plasmid_length',
        'incompatibility_group',
    ]
    with open(file_paths['mash_custom_plasmid_candidates_path'], 'w+') as candidates_file:
        writer = csv.DictWriter(candidates_file, candidates_keys, delimiter='\t', extrasaction='ignore')
        writer.writerows(mash_screen_results)


    candidates = []
    with open(file_paths['mash_custom_plasmid_candidates_path'], 'r') as candidates_file:
        reader = csv.DictReader(candidates_file, fieldnames=candidates_keys, delimiter='\t')
        for row in reader:
            candidates.append(row)

    samtools_faidx_jobs = []
    bwa_index_jobs = []
    for candidate in candidates:
        candidate_fasta_db_path = "/".join([
            mash_custom_plasmid_db,
            candidate['accession'] + ".fna"
        ])
        candidate_fasta_output_path = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".fna",
        ])
        shutil.copyfile(candidate_fasta_db_path, candidate_fasta_output_path)
        samtools_faidx_job = {
            'job_name': 'samtools_faidx',
            'native_specification': '-pe smp 2',
            'remote_command': os.path.join(job_script_path, 'samtools_faidx.sh'),
            'args': [
                "--fasta", candidate_fasta_output_path,
            ]
        }
        bwa_index_job = {
            'job_name': 'bwa_index',
            'native_specification': '-pe smp 2',
            'remote_command': os.path.join(job_script_path, 'bwa_index.sh'),
            'args': [
                "--fasta", candidate_fasta_output_path,
            ]
        }
        samtools_faidx_jobs.append(samtools_faidx_job)
        bwa_index_jobs.append(bwa_index_job)

    run_jobs(samtools_faidx_jobs + bwa_index_jobs)

    bwa_mem_jobs = []
    for candidate in candidates:
        reference_fasta = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".fna",
        ])

        bwa_mem_job = {
            'job_name': 'bwa_mem',
            'native_specification': '-pe smp 8 -shell y',
            'remote_command': os.path.join(job_script_path, 'bwa_mem.sh'),
            'args': [
                "--reference", reference_fasta,
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output", re.sub("\.fna$", "", reference_fasta) + ".sam"
            ]
        }
        bwa_mem_jobs.append(bwa_mem_job)

    run_jobs(bwa_mem_jobs)


    samtools_filter_fixmate_sort_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".sam",
        ])
        samtools_filter_fixmate_sort_job = {
            'job_name': 'samtools_filter_fixmate_sort',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_filter_fixmate_sort.sh'),
            'args': [
                "--input", alignment,
                "--flags", 1540,
                "--output", re.sub('\.sam$', '.bam', alignment),
            ]
        }
        samtools_filter_fixmate_sort_jobs.append(samtools_filter_fixmate_sort_job)

    run_jobs(samtools_filter_fixmate_sort_jobs)

    for candidate in candidates:
        sam_alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".sam",
        ])
        os.remove(sam_alignment)

    samtools_index_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".bam",
        ])
        samtools_index_job = {
            'job_name': 'samtools_index',
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_index.sh'),
            'args': [
                "--input", alignment,
            ]
        }
        samtools_index_jobs.append(samtools_index_job)

    run_jobs(samtools_index_jobs)

    samtools_depth_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".bam",
        ])
        samtools_depth_job = {
            'job_name': 'samtools_depth',
            'native_specification': '-pe smp 1',
            'remote_command': os.path.join(job_script_path, 'samtools_depth.sh'),
            'args': [
                "--input", alignment,
                "--output", re.sub('\.bam$', '.depth', alignment),
            ]
        }
        samtools_depth_jobs.append(samtools_depth_job)

    run_jobs(samtools_depth_jobs)

    for candidate in candidates:
        depth = "/".join([
            file_paths['plasmid_output_path'],
            candidate['accession'] + ".depth",
        ])
        MINIMUM_DEPTH = 10
        MINIMUM_COVERAGE_PERCENT = 95.0
        positions_above_minimum_depth = 0
        total_length = 0
        with open(depth) as depth_file:
            for line in depth_file:
                [_, position, depth] = line.split()
                total_length += 1
                if int(depth) >= MINIMUM_DEPTH:
                    positions_above_minimum_depth += 1
        percentage_above_minimum_depth = positions_above_minimum_depth / total_length
        print("percent good coverage: ", percentage_above_minimum_depth)
    
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

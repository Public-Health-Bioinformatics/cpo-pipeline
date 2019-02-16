import os
import csv
import glob
import operator
import re
import shutil
import time
from pprint import pprint

from cpo_pipeline.assembly.parsers import result_parsers
from cpo_pipeline.pipeline import run_jobs
from cpo_pipeline.plasmids import parsers

def samtools_filter_fixmate_sort_discrete_jobs(sample_id, candidates, paths):
    samtools_view_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            paths['plasmid_output_path'],
            candidate['accession'] + ".sam",
        ])
        samtools_view_job = {
            'job_name': "_".join(['samtools_view', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['samtools_view', sample_id]),
            ),
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
            paths['plasmid_output'],
            candidate['accession'] + ".mapped.dedup.bam",
        ])
        samtools_sort_job = {
            'job_name': "_".join(['samtools_sort', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['samtools_view', sample_id]),
            ),
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
            paths['plasmid_output'],
            candidate['accession'] + ".mapped.dedup.namesort.bam",
        ])
        samtools_fixmate_job = {
            'job_name': "_".join(['samtools_fixmate', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['samtools_fixmate', sample_id]),
            ),
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
            paths['plasmid_output'],
            candidate['accession'] + ".mapped.dedup.namesort.fixmate.bam",
        ])
        samtools_sort_job = {
            'job_name': "_".join(['samtools_sort', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['samtools_sort', sample_id]),
            ),
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
            paths['plasmid_output'],
            candidate['accession'] + ".mapped.dedup.namesort.fixmate.coordsort.bam",
        ])
        samtools_markdup_job = {
            'job_name': "_".join(['samtools_markdup', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['samtools_markdup', sample_id]),
            ),
            'native_specification': '-pe smp 4',
            'remote_command': os.path.join(job_script_path, 'samtools_markdup.sh'),
            'args': [
                "--input", alignment,
                "--output", re.sub("\.bam$", ".markdup.bam", alignment),
            ]
        }
        samtools_markdup_jobs.append(samtools_markdup_job)

    run_jobs(samtools_markdup_jobs)


def samtools_filter_fixmate_sort_single_job(sample_id, candidates, paths):
    samtools_filter_fixmate_sort_jobs = []
    for candidate in candidates:
        alignment = "/".join([
            paths['plasmid_output'],
            candidate['accession'] + ".sam",
        ])
        samtools_filter_fixmate_sort_job = {
            'job_name': "_".join(['samtools_filter_fixmate_sort', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['samtools_filter_fixmate_sort', sample_id]),
            ),
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
    

def custom_plasmids(sample_id, paths, queue):
    mash_jobs = [
        {
            'job_name': "_".join(['mash_screen_custom_plasmid', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['mash_screen_custom_plasmid', sample_id]),
            ),
            'native_specification': '-pe smp 8 -shell y',
            'remote_command': os.path.join(paths['job_scripts'], 'mash_screen_custom_db.sh'),
            'args': [
                "--R1", paths['reads1_fastq'],
                "--R2", paths['reads2_fastq'],
                "--min-identity", 0.996,
                "--plasmid-db-dir", os.path.join(
                    paths['mash_custom_plasmid_db'],
                    "mash",
                ),
                "--output_file", os.path.join(
                    paths['custom_plasmid_output'],
                    'mash_screen.tsv',
                )
            ],
        },
    ]
    
    run_jobs(mash_jobs)

    
    mash_screen_results = result_parsers.parse_mash_result(
        os.path.join(
            paths['custom_plasmid_output'],
            'mash_screen.tsv',
        )
    )

    custom_plasmid_db_data = {}
    for dat_file in glob.glob(os.path.join(paths['mash_custom_plasmid_db'], "data", "*.dat")):
        [dat] = parsers.custom_plasmid_db_dat_parser(dat_file)
        custom_plasmid_db_data[dat['accession']] = dat

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

    candidates_keys = [
        'identity',
        'accession',
        'circularity',
        'plasmid_length',
        'allele',
        'incompatibility_group',
    ]
    
    with open(os.path.join(paths['custom_plasmid_output'], 'candidates.tsv'), 'w+') as candidates_file:
        writer = csv.DictWriter(candidates_file, candidates_keys,
                                delimiter='\t', extrasaction='ignore')
        writer.writerows(mash_screen_results)

    candidates = []
    with open(os.path.join(paths['custom_plasmid_output'], 'candidates.tsv'), 'r') as candidates_file:
        reader = csv.DictReader(candidates_file, fieldnames=candidates_keys, delimiter='\t')
        for row in reader:
            row['fasta_path'] = os.path.join(
                paths['custom_plasmid_output'],
                'candidates',
                row['accession'] + '.fna',
            )
            candidates.append(row)

    for candidate in candidates:
        candidate_fasta_db_path = os.path.join(
            paths['mash_custom_plasmid_db'],
            candidate['accession'] + ".fna"
        )
        shutil.copyfile(candidate_fasta_db_path, candidate['fasta_path'])

    for candidate in candidates:
        candidate['database'] = 'custom'
        queue.put(candidate)
    queue.put(None)



def refseq_plasmids(sample_id, paths, queue):

    mash_jobs = [
        {
            'job_name': "_".join(['mash_screen_refseq_plasmid', sample_id]),
            'output_path': os.path.join(
                paths['plasmid_output'],
                'logs',
                "_".join(['mash_screen_refseq_plasmid', sample_id]),
            ),
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(paths['job_scripts'], 'mash_screen.sh'),
            'args': [
                "--R1", paths['reads1_fastq'],
                "--R2", paths['reads2_fastq'],
                "--queries", paths['mash_refseq_plasmid_db'],
                "--min-identity", 0.975,
                "--output_file", os.path.join(
                    paths['refseq_plasmid_output'],
                    'mash_screen.tsv',
                ),
            ],
        },
    ]
    run_jobs(mash_jobs)
    
    mash_screen_results = result_parsers.parse_mash_result(
        os.path.join(
            paths['refseq_plasmid_output'],
            'mash_screen.tsv',
        )
    )
    
    for result in mash_screen_results:
        result['accession'] = re.search('ref\|(.*)\|', result['query_id']).group(1)
        
    candidates_keys = [
        'identity',
        'accession',
    ]
    
    with open(os.path.join(paths['refseq_plasmid_output'], 'candidates.tsv'), 'w+') as candidates_file:
        writer = csv.DictWriter(candidates_file, candidates_keys,
                                delimiter='\t', extrasaction='ignore')
        writer.writerows(mash_screen_results)

    candidates = []
    with open(os.path.join(paths['refseq_plasmid_output'], 'candidates.tsv'), 'r') as candidates_file:
        reader = csv.DictReader(candidates_file, fieldnames=candidates_keys, delimiter='\t')
        for row in reader:
            row['fasta_path'] = os.path.join(
                paths['refseq_plasmid_output'],
                'candidates',
                row['accession'] + '.fna',
            )
            candidates.append(row)


    # NCBI Rate-limits downloads to 3 per second.
    # Group list of candidates into chunks of 3, then submit jobs in batches of 3 with
    # at least a 2-second delay between job submissions.
    candidates_groupby_three = [candidates[pos:pos + 3] for pos in range(0, len(candidates), 3)]
    for three_candidates in candidates_groupby_three:
        ncbi_acc_download_jobs = []
        for candidate in three_candidates:
            candidate_fasta = os.path.join(
                candidate['fasta_path']
            )

            ncbi_acc_download_job = {
                'job_name': "_".join(['ncbi_acc_download', sample_id]),
                'output_path': os.path.join(
                    paths['plasmid_output'],
                    'logs',
                    "_".join(['ncbi_acc_download', sample_id]),
                ),
                'native_specification': '-pe smp 1',
                'remote_command': os.path.join(paths['job_scripts'], 'ncbi-acc-download.sh'),
                'args': [
                    '--accession', candidate['accession'],
                    '--output_file', candidate['fasta_path']
                ]
            }
            ncbi_acc_download_jobs.append(ncbi_acc_download_job)
        time.sleep(2)
        run_jobs(ncbi_acc_download_jobs)

    for candidate in candidates:
        candidate['database'] = 'refseq'
        queue.put(candidate)
    queue.put(None)

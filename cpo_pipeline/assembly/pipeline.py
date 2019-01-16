#!/usr/bin/env python

'''
This script is a wrapper for module one of the cpo-pipeline including QC and Assembly.
It uses Mash2.0 and fastqc to check for sequence contamination, quality information and
identify a reference genome. Then assembles the reads using shovill, and assesses the
quality of the assembly with quast and busco.

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
import time
import urllib.request

import drmaa

from pkg_resources import resource_filename

from cpo_pipeline.pipeline import prepare_job
from cpo_pipeline.assembly.parsers import result_parsers

def busco_qc_check(busco_results, qc_thresholds):
    """
    BUSCO PASS CRITERIA:
    1. complete singles > 90% of total genes
    2. complte duplicates < 90% of total genes
    Args:
        busco_results (dict): Busco results
        qc_thresholds (dict): Threshold values for determining QC pass/fail
    Returns:
        boolean: Assembly passes our BUSCO quality criteria
    """
    complete_single = busco_results['complete_single']
    complete_duplicate = busco_results['complete_duplicate']
    total = busco_results['total']
    busco_complete_single_cutoff = qc_thresholds['busco_complete_single_cutoff']
    busco_complete_duplicate_cutoff = qc_thresholds['busco_complete_duplicate_cutoff']

    return bool((complete_single / total) >= busco_complete_single_cutoff and \
                (complete_duplicate / total) <= busco_complete_duplicate_cutoff)


def quast_qc_check(quast_results, qc_thresholds):
    """
    QUAST PASS CRITERIA:
    1. total length vs reference length +-10%
    2. percent gc versus reference percent gc +- 5%
    3. genome fraction percent > 90
    Args:
        quast_results (dict): Quast results
        qc_thresholds (dict): Threshold values for determining QC pass/fail
    Returns:
        boolean: Assembly passes our QUAST quality criteria
    """
    total_length = quast_results['total_length']
    reference_length = quast_results['reference_length']
    assembly_percent_gc = quast_results['percent_GC']
    reference_percent_gc = quast_results['reference_percent_GC']
    genome_fraction_percent = quast_results['genome_fraction_percent']
    percent_gc_cutoff = qc_thresholds['quast_percent_gc_cutoff']
    assembly_length_cutoff = qc_thresholds['quast_assembly_length_cutoff']
    genome_fraction_percent_cutoff = qc_thresholds['quast_genome_fraction_percent_cutoff']
    return bool(total_length <= (reference_length * (1 + assembly_length_cutoff)) and \
                total_length >= (reference_length * (1 - assembly_length_cutoff)) and \
                assembly_percent_gc <= (reference_percent_gc * (1 + percent_gc_cutoff)) and \
                assembly_percent_gc >= (reference_percent_gc * (1 - percent_gc_cutoff)) and \
                genome_fraction_percent >= genome_fraction_percent_cutoff)


def fastqc_qc_check(fastqc_results):
    """
    Args:
        fastqc_results (dict): FastQC results
    Returns:
        boolean: Sequence data passes our FastQC quality criteria
    """
    return bool(fastqc_results["basic_statistics"] == "PASS" and \
                fastqc_results["per_base_sequence_quality"] == "PASS" and \
                fastqc_results["sequence_length_distribution"] == "PASS")


def download_mash_hit(mash_hit, download_path):
    """
    Given a mash_hit, download the query sequence from NCBI FTP servers
    Will fail if the download_path doesn't exist.
    Args:
        mash_hit(dict):
        download_path(str):
    Returns:
        (void)
    """

    def mash_query_id_to_ncbi_ftp_path(query_id):
        """
        Args:
            query_id (str): Mash query ID (column 5 of mash screen report)
        Returns:
            list: Directory names used to locate reference genome
                  on ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
        For example:
            "GCF/001/022/155"
        """
        prefix = query_id.split('_')[0]
        digits = query_id.split('_')[1].split('.')[0]
        path_list = [prefix] + [digits[i:i+3] for i in range(0, len(digits), 3)]

        return "/".join(path_list)

    query_id = mash_hit['query_id']
    ncbi_ftp_path = mash_query_id_to_ncbi_ftp_path(query_id)
    assembly = query_id[:query_id.find("_genomic.fna.gz")]

    ncbi_ftp_server_base = "ftp://ftp.ncbi.nlm.nih.gov"
    fasta_url = "/".join([
        ncbi_ftp_server_base, "genomes", "all",
        ncbi_ftp_path,
        assembly,
        query_id
    ])
    assembly_stat_url = "/".join([
        ncbi_ftp_server_base, "genomes", "all",
        ncbi_ftp_path,
        assembly,
        assembly + "_assembly_stats.txt"
    ])

    #fetch the files
    urllib.request.urlretrieve(fasta_url, "/".join([download_path, query_id]))
    urllib.request.urlretrieve(assembly_stat_url,
                               "/".join([download_path, assembly + "_assembly_stats.txt"]))

def prepare_output_directories(output_dir, sample_id):
    """
    Prepare output sub-directories before running analyses
    Args:
        output_dir(str): Base output directory
        sample_id(str): Sample ID
    Returns:
        (void)
    """
    output_subdirs = [
        "/".join([output_dir, sample_id, subdir])
        for subdir in [
            'pre-assembly_qc',
            'assembly',
            'post-assembly_qc',
            'reference'
        ]
    ]
    for output_subdir in output_subdirs:
        try:
            os.makedirs(output_subdir)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

def main(parser, config):
    """
    main entrypoint
    Args:
        parser():
        config():
    Returns:
        (void)
    """
    if not parser:
        script_name = os.path.basename(os.path.realpath(sys.argv[0]))
        parser = argparse.ArgumentParser(prog=script_name)
        parser.add_argument("-i", "--ID", dest="sample_id",
                            help="identifier of the isolate")
        parser.add_argument("-1", "--R1", dest="reads1_fastq",
                            help="absolute file path forward read (R1)")
        parser.add_argument("-2", "--R2", dest="reads2_fastq",
                            help="absolute file path to reverse read (R2)")
        parser.add_argument("-o", "--output", dest="output", default='./',
                            help="absolute path to output folder")
    parser.add_argument("-g", "--mash-genomedb", dest="mash_genome_db",
                        help="absolute path to mash reference database")
    parser.add_argument("-p", "--mash-plasmiddb", dest="mash_plasmid_db",
                        help="absolute path to mash reference database")
    parser.add_argument("-b", "--busco-db", dest="busco_db",
                        help="absolute path to busco reference database")

    args = parser.parse_args()

    if not config:
        config = configparser.ConfigParser()
        config_file = resource_filename('data', 'config.ini')
        config.read(config_file)
    
    if args.mash_genome_db and not config['databases']['mash_genome_db']:
        mash_genome_db = args.mash_genome_db
    else:
        mash_genome_db = config['databases']['mash_genome_db']
    if args.mash_plasmid_db and not config['databases']['mash_plasmid_db']:
        mash_plasmid_db = args.mash_plasmid_db
    else:
        mash_plasmid_db = config['databases']['mash_plasmid_db']
    if args.busco_db and not config['databases']['busco_db']:
        busco_db = args.busco_db
    else:
        busco_db = config['databases']['busco_db']
    sample_id = args.sample_id
    reads1_fastq = args.reads1_fastq
    reads2_fastq = args.reads2_fastq
    output_dir = args.output
    
    prepare_output_directories(output_dir, sample_id)

    #dictionary to store QC PASS/FAIL flags
    qc_verdicts = {
        "multiple_species_contamination": None,
        "fastq_contains_plasmids": None,
        "acceptable_coverage": None,
        "acceptable_fastqc_forward": None,
        "acceptable_fastqc_reverse": None,
        "acceptable_quast_assembly_metrics": None,
        "acceptable_busco_assembly_metrics": None
    }

    qc_thresholds = {
        # genome mash will include all hits with scores (top hit score - $thisvalue)
        "mash_hits_genome_score_cutoff": 300,
        # plasmid mash will include all hits with scores (top hit score - $thisvalue)
        "mash_hits_plasmid_score_cutoff": 100,
        # sequencing coverage greater than ($thisvalue) will pass the QC
        "coverage_cutoff": 30,
        # QUAST QC: assembly length within +-($thisvalue) percent
        # in reference to reference length will pass the QC
        "quast_assembly_length_cutoff": 0.10,
        # QUAST QC: percent GC within +-($thisvalue) percent in reference
        # to reference percent GC will pass the QC
        "quast_percent_gc_cutoff":0.05,
        # QUAST QC: genome_fraction_percent greater than ($thisvalue) will pass the QC
        "quast_genome_fraction_percent_cutoff":0.90,
        # BUSCO QC: complete single genes greater than ($thisvalue) percent will pass the QC
        "busco_complete_single_cutoff":0.90,
        # BUSCO QC: complete duplicate genes less than ($thisvalue) percent will pass the QC
        "busco_complete_duplicate_cutoff":0.10
    }

    file_paths = {
        "mash_genome_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "mashscreen.genome.tsv"]),
        "mash_plasmid_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "mashscreen.plasmid.tsv"]),
        "fastqc_output_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "fastqc"]),
        "totalbp_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "totalbp"]),
        "reference_genome_path": "/".join([output_dir, sample_id, "reference"]),
        "assembly_path": "/".join([output_dir, sample_id, "assembly"]),
        "busco_path": "/".join([output_dir, sample_id, "post-assembly_qc", "busco"]),
        "quast_path": "/".join([output_dir, sample_id, "post-assembly_qc", "quast"]),
    }

    print(str(datetime.datetime.now()))
    print("sample_id: " + sample_id)
    print("reads1_fastq: " + reads1_fastq)
    print("reads2_fastq: " + reads2_fastq)

    print("step 1: preassembly QC")

    job_script_path = resource_filename('data', 'job_scripts')

    pre_assembly_qc_jobs = [
        {
            'job_name': 'mash_screen',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_screen.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--queries", mash_genome_db,
                "--output_file", file_paths['mash_genome_path']
            ],
        },
        {
            'job_name': 'mash_screen',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_screen.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--queries", mash_plasmid_db,
                "--output_file", file_paths['mash_plasmid_path']
            ],
        },
        {
            'job_name': 'fastqc',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'fastqc.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output_dir", file_paths['fastqc_output_path']
            ],
        },
        {
            'job_name': 'seqtk_totalbp',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'seqtk_totalbp.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output_file", file_paths['totalbp_path']
            ],
        }
    ]

    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in pre_assembly_qc_jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    print("Parsing the QC results")
    #parse genome mash results
    mash_hits = result_parsers.parse_mash_result(file_paths["mash_genome_path"])
    mash_hits = sorted(mash_hits, key=lambda k: k['identity'], reverse=True)
    # 'shared_hashes' field is string in format '935/1000'
    # Threshold is 300 below highest numerator (ie. '935/100' -> 635)
    mash_hits_score_threshold = int(mash_hits[0]['shared_hashes'].split("/")[0]) - \
        int(qc_thresholds["mash_hits_genome_score_cutoff"])

    def score_above_threshold(mash_result, score_threshold):
        score = int(mash_result['shared_hashes'].split("/")[0])
        return bool(score >= score_threshold and \
                    mash_result['query_comment'].find("phiX") == -1)

    filtered_mash_hits = list(filter(
        lambda x: score_above_threshold(x, mash_hits_score_threshold),
        mash_hits))

    # parse plasmid mash
    mash_plasmid_hits = result_parsers.parse_mash_result(file_paths["mash_plasmid_path"])
    mash_plasmid_hits = sorted(mash_plasmid_hits, key=lambda k: k['identity'], reverse=True)
    # 'shared_hashes' field is string in format '935/1000'
    # Threshold is 100 below highest numerator (ie. '935/100' -> 835)
    mash_plasmid_hits_score_threshold = int(mash_plasmid_hits[0]['shared_hashes'].split("/")[0]) - \
        int(qc_thresholds["mash_hits_plasmid_score_cutoff"])
    filtered_mash_plasmid_hits = list(filter(
        lambda x: score_above_threshold(x, mash_plasmid_hits_score_threshold),
        mash_plasmid_hits))

    # parse fastqc
    fastqc_r1_result = result_parsers.parse_fastqc_result(
        file_paths['fastqc_output_path'] + "/" + sample_id + "_R1_fastqc/summary.txt"
    )
    fastqc_r2_result = result_parsers.parse_fastqc_result(
        file_paths['fastqc_output_path'] + "/" + sample_id + "_R2_fastqc/summary.txt"
    )

    #all the qC result are parsed now, lets do some QC logic
    #look at mash results first
    if len(filtered_mash_hits) > 1:
        qc_verdicts["multiple_species_contamination"] = True
    else:
        qc_verdicts["multiple_species_contamination"] = False


    if filtered_mash_plasmid_hits:
        qc_verdicts["fastq_contains_plasmids"] = True
    else:
        qc_verdicts["fastq_contains_plasmids"] = False

    #look at fastqc results
    qc_verdicts["acceptable_fastqc_forward"] = fastqc_qc_check(fastqc_r1_result)
    qc_verdicts["acceptable_fastqc_reverse"] = fastqc_qc_check(fastqc_r2_result)

    #download a reference genome
    print("Downloading reference genomes")

    reference_genomes = []
    # build the save paths
    try:
        os.makedirs(file_paths['reference_genome_path'])
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
    if not qc_verdicts["multiple_species_contamination"]:
        for mash_hit in filtered_mash_hits:
            if not re.match("phiX", mash_hit['query_comment']):
                download_mash_hit(mash_hit, file_paths['reference_genome_path'])
                reference_genomes.append(mash_hit['query_id'])

    else: #throw an error if it contains contaminations
        print("Contaminated Genome assembly...resequencing required")
        raise Exception("contamination and mislabeling...crashing")

    #check to make sure we ONLY have ONE reference.
    if len(reference_genomes) > 1:
        raise Exception("there are multiple reference genomes")
    elif not reference_genomes:
        raise Exception("no reference genome identified")

    # now we estimate our coverage using total reads and expected genome size
    expected_genome_size = result_parsers.parse_reference_genome_stats(
        glob.glob(file_paths["reference_genome_path"] + "/*_assembly_stats.txt")[0]
    )
    total_bp = result_parsers.parse_total_bp(file_paths["totalbp_path"])

    coverage = total_bp / expected_genome_size

    if coverage >= int(qc_thresholds["coverage_cutoff"]):
        qc_verdicts["acceptable_coverage"] = True


    print("step 2: genome assembly and QC")

    assembly_jobs = [
        {
            'job_name': 'shovill',
            'native_specification': '-pe smp 16 -l h_vmem=4G',
            'remote_command': os.path.join(job_script_path, 'shovill.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--mincov", "3",
                "--minlen", "500",
                "--output_dir", file_paths['assembly_path']
            ],
        }
    ]

    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in assembly_jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    post_assembly_qc_jobs = [
        {
            'job_name': 'busco',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'busco.sh'),
            'args': [
                "--input", "/".join([file_paths['assembly_path'], "contigs.fa"]),
                "--database", busco_db,
                "--output_dir", file_paths['busco_path']
            ]
        },
        {
            'job_name': 'quast',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'quast.sh'),
            'args': [
                "--input", "/".join([file_paths['assembly_path'], "contigs.fa"]),
                "--reference_genome", "/".join(
                    [file_paths['reference_genome_path'], reference_genomes[0]]
                ),
                "--output_dir", file_paths['quast_path']
            ]
        },
    ]

    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in post_assembly_qc_jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    print("Parsing assembly results")
    busco_results = result_parsers.parse_busco_result(
        file_paths["busco_path"] + "/run_busco/short_summary_busco.txt"
    )
    quast_results = result_parsers.parse_quast_result(
        file_paths["quast_path"] + "/report.txt"
    )

    qc_verdicts["acceptable_busco_assembly_metrics"] = busco_qc_check(busco_results, qc_thresholds)
    qc_verdicts["acceptable_quast_assembly_metrics"] = quast_qc_check(quast_results, qc_thresholds)

    #print QC results to screen
    print("total bases: " + str(total_bp))
    print("expected genome size: " + str(expected_genome_size))
    print("coverage: " + str(coverage))
    print("")
    for key, value in qc_verdicts.items():
        print(str(key) + ": " + str(value))

    return "/".join([file_paths['assembly_path'], "contigs.fa"])

if __name__ == "__main__":
    START = time.time()
    print("Starting workflow...")
    main(None, None)
    END = time.time()
    print("Finished!\nThe analysis used: " + str(END - START) + " seconds")

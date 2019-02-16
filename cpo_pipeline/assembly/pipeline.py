#!/usr/bin/env python

'''
This script is a wrapper for module one of the cpo-pipeline including QC and Assembly.
It uses Mash2.0 and fastqc to check for sequence contamination, quality information and
identify a reference genome. Then assembles the reads using shovill, and assesses the
quality of the assembly with quast.

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

from cpo_pipeline.pipeline import prepare_job, run_jobs
from cpo_pipeline.assembly.parsers import result_parsers
from cpo_pipeline.assembly.parsers import input_parsers

# If the user provides an unrecognized NCBI taxid
# then use this as the estimated genome size
DEFAULT_ESTIMATED_GENOME_SIZE = 5000000

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


def quast_qc_check(quast_results, estimated_genome_size):
    """
    QUAST PASS CRITERIA:
    1. total length within +/- 10% of expected genome size
    Args:
        quast_results (dict): Quast results
        qc_thresholds (dict): Threshold values for determining QC pass/fail
    Returns:
        boolean: Assembly passes our QUAST quality criteria
    """
    total_length = quast_results['total_length']
    return bool(total_length <= (estimated_genome_size * 1.1) and \
                total_length >= (estimated_genome_size * 0.9))


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
            'reference',
        ]
    ]
    for output_subdir in output_subdirs:
        try:
            os.makedirs(output_subdir)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

def get_estimated_genome_size(estimated_genome_sizes, ncbi_taxonomy_id):
    """
    Args:
        estimated_genome_sizes(dict): parsed estimated_genome_sizes.tsv file.
        ncbi_taxonomy_id(str): NCBI Taxonomy ID
    Returns:
        (int) Estimated genome size
    """
    try:
        [estimated_genome_size] = [x['estimated_genome_size'] for x in estimated_genome_sizes if
                                   x['ncbi_taxonomy_id'] == ncbi_taxonomy_id]
    except ValueError:
        estimated_genome_size = DEFAULT_ESTIMATED_GENOME_SIZE

    return estimated_genome_size


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

    if args.mash_genome_db and not config['databases']['mash_genome_db']:
        mash_genome_db = args.mash_genome_db
    else:
        mash_genome_db = config['databases']['mash_genome_db']

    sample_id = args.sample_id
    reads1_fastq = args.reads1_fastq
    reads2_fastq = args.reads2_fastq
    output_dir = args.outdir

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
        # BUSCO QC: complete single genes greater than ($thisvalue) percent will pass the QC
        "busco_complete_single_cutoff":0.90,
        # BUSCO QC: complete duplicate genes less than ($thisvalue) percent will pass the QC
        "busco_complete_duplicate_cutoff":0.10
    }

    file_paths = {
        "mash_genome_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "mashscreen.genome.tsv"]),
        "fastqc_output_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "fastqc"]),
        "totalbp_path": "/".join([output_dir, sample_id, "pre-assembly_qc", "totalbp"]),
        "reference_genome_path": "/".join([output_dir, sample_id, "reference"]),
        "assembly_path": "/".join([output_dir, sample_id, "assembly"]),
        "quast_path": "/".join([output_dir, sample_id, "post-assembly_qc", "quast"]),
    }

    print(str(datetime.datetime.now()))
    print("sample_id: " + sample_id)
    print("reads1_fastq: " + reads1_fastq)
    print("reads2_fastq: " + reads2_fastq)

    print("step 1: preassembly QC")

    job_script_path = resource_filename('data', 'job_scripts')
    estimated_genome_sizes_path = resource_filename('data', 'estimated_genome_sizes.tsv')
    estimated_genome_sizes = input_parsers.parse_estimated_genome_sizes(estimated_genome_sizes_path)


    pre_assembly_qc_jobs = [
        {
            'job_name': "_".join(['mash_screen', sample_id]),
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
            'job_name': "_".join(['fastqc', sample_id]),
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'fastqc.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output_dir", file_paths['fastqc_output_path']
            ],
        },
        {
            'job_name': "_".join(['seqtk_totalbp', sample_id]),
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'seqtk_totalbp.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output_file", file_paths['totalbp_path']
            ],
        }
    ]

    run_jobs(pre_assembly_qc_jobs)

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

    # parse fastqc
    fastqc_r1_result = result_parsers.parse_fastqc_result(
        glob.glob(
            "/".join([
                file_paths['fastqc_output_path'],
                sample_id + "*_R1_*" + "fastqc",
                "summary.txt"
            ])
        )[0]
    )
    fastqc_r2_result = result_parsers.parse_fastqc_result(
        glob.glob(
            "/".join([
                file_paths['fastqc_output_path'],
                sample_id + "*_R2_*" + "fastqc",
                "summary.txt"
            ])
        )[0]
    )

    #all the qC result are parsed now, lets do some QC logic
    #look at mash results first
    if len(filtered_mash_hits) > 1:
        qc_verdicts["multiple_species_contamination"] = True
    else:
        qc_verdicts["multiple_species_contamination"] = False


    #look at fastqc results
    qc_verdicts["acceptable_fastqc_forward"] = fastqc_qc_check(fastqc_r1_result)
    qc_verdicts["acceptable_fastqc_reverse"] = fastqc_qc_check(fastqc_r2_result)


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

    # If the user passes an expected organism NCBI taxonomy ID, then
    # use that to estimate the genome size. Otherwise, use the downloaded reference.
    if args.expected_organism_ncbi_taxid:
        estimated_genome_size = get_estimated_genome_size(estimated_genome_sizes, args.expected_organism_ncbi_taxid)
    else:
        estimated_genome_size = result_parsers.parse_reference_genome_stats(
            glob.glob(file_paths["reference_genome_path"] + "/*_assembly_stats.txt")[0]
        )

    total_bp = result_parsers.parse_total_bp(file_paths["totalbp_path"])

    coverage = total_bp / estimated_genome_size

    if coverage >= int(qc_thresholds["coverage_cutoff"]):
        qc_verdicts["acceptable_coverage"] = True


    print("step 2: genome assembly and QC")

    assembly_jobs = [
        {
            'job_name': "_".join(['shovill', sample_id]),
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

    run_jobs(assembly_jobs)

    post_assembly_qc_jobs = [
        {
            'job_name': "_".join(['quast', sample_id]),
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'quast.sh'),
            'args': [
                "--input", "/".join([file_paths['assembly_path'], "contigs.fa"]),
                "--outdir", file_paths['quast_path']
            ]
        },
    ]

    run_jobs(post_assembly_qc_jobs)

    busco_results = result_parsers.parse_busco_result(
        file_paths["quast_path"] + "/busco_stats/short_summary_contigs.txt"
    )
    quast_results = result_parsers.parse_quast_result(
        file_paths["quast_path"] + "/report.txt"
    )

    qc_verdicts["acceptable_busco_assembly_metrics"] = busco_qc_check(busco_results, qc_thresholds)
    qc_verdicts["acceptable_quast_assembly_metrics"] = quast_qc_check(quast_results, estimated_genome_size)

    return "/".join([file_paths['assembly_path'], "contigs.fa"])

if __name__ == "__main__":
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name, description='')
    parser.add_argument("-i", "--ID", dest="sample_id",
                        help="identifier of the isolate", required=True)
    parser.add_argument("-e", "--expected-organism-taxid", dest="expected_organism_ncbi_taxid",
                        help="Expected organism NCBI Taxonomy ID")
    parser.add_argument("-1", "--R1", dest="reads1_fastq",
                        help="absolute file path forward read (R1)", required=True)
    parser.add_argument("-2", "--R2", dest="reads2_fastq",
                        help="absolute file path to reverse read (R2)", required=True)
    parser.add_argument("-o", "--outdir", dest="output", default='./',
                        help="absolute path to output folder")
    parser.add_argument('-c', '--config', dest='config_file',
                        default=resource_filename('data', 'config.ini'),
                        help='Config File', required=False)
    args = parser.parse_args()
    main(args)

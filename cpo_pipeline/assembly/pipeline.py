#!/usr/bin/env python

import argparse
import configparser
import csv
import datetime
import drmaa
import errno
import glob
import logging
import os
import re
import structlog
import sys
import urllib.request
import uuid

from pkg_resources import resource_filename

import cpo_pipeline
from cpo_pipeline.drmaa import prepare_job, run_jobs
from cpo_pipeline.logging import now
from cpo_pipeline.assembly.parsers import result_parsers
from cpo_pipeline.assembly.parsers import input_parsers
from cpo_pipeline.assembly import quality_control as qc


# If the user provides an unrecognized NCBI taxid
# then use this as the estimated genome size
DEFAULT_ESTIMATED_GENOME_SIZE = 5000000

logger = structlog.get_logger()

def download_refseq_reference(reference_id, download_path):
    """
    Given a mash_hit, download the query sequence from NCBI FTP servers
    Will fail if the download_path doesn't exist.
    Args:
        mash_hit(dict):
        download_path(str):
    Returns:
        (void)
    """

    def mash_reference_id_to_ncbi_ftp_path(reference_id):
        """
        Args:
            query_id (str): Mash reference ID (column 1 of mash dist report)
        Returns:
            list: Directory names used to locate reference genome
                  on ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
        For example:
            "GCF/001/022/155"
        """
        prefix = reference_id.split('_')[0]
        digits = reference_id.split('_')[1].split('.')[0]
        path_list = [prefix] + [digits[i:i+3] for i in range(0, len(digits), 3)]

        return "/".join(path_list)

    ncbi_ftp_path = mash_reference_id_to_ncbi_ftp_path(reference_id)
    assembly = reference_id[:reference_id.find("_genomic.fna.gz")]

    ncbi_ftp_server_base = "ftp://ftp.ncbi.nlm.nih.gov"
    fasta_url = "/".join([
        ncbi_ftp_server_base, "genomes", "all",
        ncbi_ftp_path,
        assembly,
        reference_id
    ])
    assembly_stat_url = "/".join([
        ncbi_ftp_server_base, "genomes", "all",
        ncbi_ftp_path,
        assembly,
        assembly + "_assembly_stats.txt"
    ])

    #fetch the files
    try:
        urllib.request.urlretrieve(fasta_url, "/".join([download_path, reference_id]))
        logger.info(
            "file_downloaded",
            timestamp=str(now()),
            url=fasta_url,
        )
    except Exception as e:
        logging.error(
            "download_failed",
            timestamp=str(now()),
            url=fasta_url,
        )
    try:
        urllib.request.urlretrieve(assembly_stat_url,
                               "/".join([download_path, assembly + "_assembly_stats.txt"]))
        logger.info(
            "file_downloaded",
            timestamp=str(now()),
            url=assembly_stat_url,
        )
    except Exception as e:
        logging.error(
            "download_failed",
            timestamp=str(now()),
            url=assembly_stat_url,
        )

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
        os.path.join(output_dir, sample_id, subdir)
        for subdir in [
            'logs',
            'pre-assembly_qc',
            'assembly',
            'post-assembly_qc',
            'reference',
        ]
    ]
    for output_subdir in output_subdirs:
        try:
            os.makedirs(output_subdir)
        except OSError as e:
            if e.errno != errno.EEXIST:
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

    try:
        mash_genome_db = args.mash_genome_db
    except AttributeError:
        try:
            mash_genome_db = config['databases']['mash_genome_db']
            if not os.path.exists(mash_genome_db):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), mash_genome_db
                )
            logger.info(
                "configuration_loaded",
                timestamp=str(now()),
                configuration_attribute="databases/mash_genome_db",
                configuration_value=mash_genome_db,
            )
        except Exception as e:
            logger.error(
                "configuration_failed",
                timestamp=str(now()),
                configuration_attribute="databases/mash_genome_db",
                error_message=str(e),
            )

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

    paths = {
        "output_dir": output_dir,
        'logs': os.path.join(
            output_dir,
            sample_id,
            'logs',
        ),
        "mash_genome_path": os.path.join(
            output_dir,
            sample_id,
            "pre-assembly_qc",
            "mash_dist.genome.tsv"
        ),
        "fastqc_output_path": os.path.join(
            output_dir,
            sample_id,
            "pre-assembly_qc",
            "fastqc"
        ),
        "totalbp_path": os.path.join(
            output_dir,
            sample_id,
            "pre-assembly_qc",
            "totalbp"
        ),
        "estimated_coverage_stats_path": os.path.join(
            output_dir,
            sample_id,
            "pre-assembly_qc",
            "estimated_coverage_stats.tsv"
        ),
        "reference_genome_path": os.path.join(
            output_dir,
            sample_id,
            "reference"
        ),
        "assembly_output": os.path.join(
            output_dir,
            sample_id,
            "assembly"
        ),
        "quast_path": os.path.join(
            output_dir,
            sample_id,
            "post-assembly_qc",
            "quast"
        ),
    }


    job_script_path = resource_filename('data', 'job_scripts')
    estimated_genome_sizes_path = resource_filename('data', 'estimated_genome_sizes.tsv')
    estimated_genome_sizes = input_parsers.parse_estimated_genome_sizes(estimated_genome_sizes_path)


    pre_assembly_qc_jobs = [
        {
            'job_name': "_".join(['mash_dist_sort_head', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_dist_sort_head.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--queries", mash_genome_db,
                "--output_file", paths['mash_genome_path']
            ],
        },
        {
            'job_name': "_".join(['fastqc', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'fastqc.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output_dir", paths['fastqc_output_path']
            ],
        },
        {
            'job_name': "_".join(['seqtk_totalbp', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'seqtk_totalbp.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--output_file", paths['totalbp_path']
            ],
        }
    ]

    run_jobs(pre_assembly_qc_jobs)

    #parse genome mash results
    mash_dist_results = []
    try:
        mash_dist_results = result_parsers.parse_mash_dist_result(paths["mash_genome_path"])
        logger.info(
            "parsed_result_file",
            timestamp=str(now()),
            filename=os.path.abspath(paths["mash_genome_path"]),
            closest_match_reference_id=mash_dist_results[0]['reference_id'],
        )
    except Exception as e:
        logger.info(
            "result_parsing_failed",
            timestamp=str(now()),
            filename=os.path.abspath(paths["mash_genome_path"]),
            error_message=e.message,
        )


    # parse fastqc
    fastqc_results = {}
    for read in ["R1", "R2"]:
        try:
            [fastqc_result_summary_path] = glob.glob(
                os.path.join(
                    paths['fastqc_output_path'],
                    "*_" + read + "_*" + "fastqc",
                    'summary.txt'
                )
            )
            fastqc_results[read] = result_parsers.parse_fastqc_result(
                fastqc_result_summary_path
            )
            logger.info(
                "parsed_result_file",
                timestamp=str(now()),
                filename=os.path.abspath(fastqc_result_summary_path),
                summary=fastqc_results[read],
            )
        except Exception as e:
            logger.error(
                "result_parsing_failed",
                timestamp=str(now()),
                filename=fastqc_result_summary_path
            )
            fastqc_results["R1"] = {
                "basic_statistics": "FAILED_TO_PARSE",
                "per_base_sequence_quality": "FAILED_TO_PARSE",
                "per_tile_sequence_quality": "FAILED_TO_PARSE",
                "per_sequence_quality_scores": "FAILED_TO_PARSE",
                "per_base_sequence_content": "FAILED_TO_PARSE",
                "per_sequence_gc_content": "FAILED_TO_PARSE",
                "per_base_n_content": "FAILED_TO_PARSE",
                "sequence_length_distribution": "FAILED_TO_PARSE",
                "sequence_duplication_levels": "FAILED_TO_PARSE",
                "overrepresented_sequences": "FAILED_TO_PARSE",
                "adapter_content": "FAILED_TO_PARSE",
            }

            fastqc_results["R2"] = {
                "basic_statistics": "FAILED_TO_PARSE",
                "per_base_sequence_quality": "FAILED_TO_PARSE",
                "per_tile_sequence_quality": "FAILED_TO_PARSE",
                "per_sequence_quality_scores": "FAILED_TO_PARSE",
                "per_base_sequence_content": "FAILED_TO_PARSE",
                "per_sequence_gc_content": "FAILED_TO_PARSE",
                "per_base_n_content": "FAILED_TO_PARSE",
                "sequence_length_distribution": "FAILED_TO_PARSE",
                "sequence_duplication_levels": "FAILED_TO_PARSE",
                "overrepresented_sequences": "FAILED_TO_PARSE",
                "adapter_content": "FAILED_TO_PARSE",
            }


    #look at fastqc results
    qc_verdicts["acceptable_fastqc_forward"] = qc.fastqc_qc_check(fastqc_results["R1"])
    qc_verdicts["acceptable_fastqc_reverse"] = qc.fastqc_qc_check(fastqc_results["R2"])

    try:
        reference_genome = mash_dist_results[0]['reference_id']
    except Exception as e:
        logger.error(
            "failed_quality_control_check",
            timestamp=str(now()),
            qc_check_failed="determine_reference_sequence",
            error_message=e.message,
        )
        
    # build the save paths
    try:
        os.makedirs(paths['reference_genome_path'])
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise

    download_refseq_reference(reference_genome, paths['reference_genome_path'])


    # If the user passes an expected organism NCBI taxonomy ID, then
    # use that to estimate the genome size. Otherwise, use the downloaded reference.
    estimated_genome_size = DEFAULT_ESTIMATED_GENOME_SIZE
    if args.expected_organism_ncbi_taxid:
        estimated_genome_size = get_estimated_genome_size(estimated_genome_sizes, args.expected_organism_ncbi_taxid)
    else:
        try:
            [reference_genome_assembly_stats_path] = glob.glob(
                paths["reference_genome_path"] + "/*_assembly_stats.txt"
            )
        except ValueError:
            logger.error(
                "result_parsing_failed",
                timestamp=str(now()),
                filename=str(os.path.abspath(paths["reference_genome_path"])) + "/*_assembly_stats.txt",
            )

        try:
            reference_genome_assembly_stats = result_parsers.parse_reference_genome_assembly_stats(
                reference_genome_assembly_stats_path
            )
            logger.info(
                "parsed_result_file",
                timestamp=str(now()),
                filename=os.path.abspath(reference_genome_assembly_stats_path),
                total_length=reference_genome_assembly_stats['total_length'],
                contig_count=reference_genome_assembly_stats['contig_count'],
                contig_N50=reference_genome_assembly_stats['contig_N50'],
                organism_name=reference_genome_assembly_stats['organism_name'],
                infraspecific_name=reference_genome_assembly_stats['infraspecific_name'],
                ncbi_taxonomy_id=reference_genome_assembly_stats['taxid'],
                refseq_assembly_accession=reference_genome_assembly_stats['refseq_assembly_accession'],
            )
            estimated_genome_size = reference_genome_assembly_stats['total_length']
        except Exception as e:
            logger.error(
                "result_parsing_failed",
                timestamp=str(now()),
                filename=os.path.abspath(reference_genome_assembly_stats_path),
                error_message=e.message,
            )

    total_bp = result_parsers.parse_total_bp(paths["totalbp_path"])
    logger.info(
        "parsed_result_file",
        timestamp=str(now()),
        filename=os.path.abspath(paths["totalbp_path"]),
        total_bp=total_bp,
    )

    estimated_depth_of_coverage = total_bp / estimated_genome_size


    if estimated_depth_of_coverage >= int(qc_thresholds["coverage_cutoff"]):
        qc_verdicts["acceptable_coverage"] = True

    estimated_coverage_stats_headers = [
        'sample_id',
        'total_bp',
        'estimated_genome_size',
        'estimated_depth_of_coverage',
    ]

    with open(paths['estimated_coverage_stats_path'], 'w+') as f:
        writer = csv.DictWriter(f, fieldnames=estimated_coverage_stats_headers, delimiter='\t')
        writer.writeheader()
        writer.writerow({
            'sample_id': sample_id,
            'total_bp': int(total_bp),
            'estimated_genome_size': int(estimated_genome_size),
            'estimated_depth_of_coverage': round(estimated_depth_of_coverage, 4),
        })

    assembly_jobs = [
        {
            'job_name': "_".join(['shovill', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 16 -l h_vmem=4G',
            'remote_command': os.path.join(job_script_path, 'shovill.sh'),
            'args': [
                "--R1", reads1_fastq,
                "--R2", reads2_fastq,
                "--mincov", "3",
                "--minlen", "500",
                "--output_dir", paths['assembly_output']
            ],
        }
    ]

    run_jobs(assembly_jobs)

    post_assembly_qc_jobs = [
        {
            'job_name': "_".join(['quast', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'quast.sh'),
            'args': [
                "--input", os.path.join(paths['assembly_output'], "contigs.fa"),
                "--outdir", paths['quast_path']
            ]
        },
    ]

    run_jobs(post_assembly_qc_jobs)

    busco_short_summary_contigs_path = os.path.abspath(
        paths["quast_path"] + "/busco_stats/short_summary_contigs.txt"
    )
    busco_results = result_parsers.parse_busco_result(
        busco_short_summary_contigs_path
    )
    logger.info(
        "parsed_result_file",
        timestamp=str(now()),
        filename=os.path.abspath(busco_short_summary_contigs_path),
        busco_results=busco_results
    )
    quast_report_path = os.path.abspath(paths["quast_path"] + "/report.txt")
    quast_results = result_parsers.parse_quast_result(
        quast_report_path
    )
    logger.info(
        "parsed_result_file",
        timestamp=str(now()),
        filename=os.path.abspath(quast_report_path),
        num_contigs=quast_results["num_contigs"],
        N50=quast_results["N50"],
    )

    qc_verdicts["acceptable_busco_assembly_metrics"] = qc.busco_qc_check(busco_results, qc_thresholds)
    qc_verdicts["acceptable_quast_assembly_metrics"] = qc.quast_qc_check(quast_results, estimated_genome_size)


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
    parser.add_argument("-o", "--outdir", dest="outdir", default='./',
                        help="absolute path to output folder")
    parser.add_argument('-c', '--config', dest='config_file',
                        default=resource_filename('data', 'config.ini'),
                        help='Config File', required=False)

    args = parser.parse_args()
    
    logging.basicConfig(
        format="%(message)s",
        stream=sys.stdout,
        level=logging.DEBUG,
    )

    structlog.configure_once(
        processors=[
            structlog.stdlib.add_log_level,
            structlog.processors.JSONRenderer()
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        context_class=structlog.threadlocal.wrap_dict(dict),
    )

    logger = structlog.get_logger(
        analysis_id=str(uuid.uuid4()),
        sample_id=args.sample_id,
        pipeline_version=cpo_pipeline.__version__,
    )
    
    logger.info(
        "analysis_started",
        timestamp=str(now()),
    )
    
    main(args)

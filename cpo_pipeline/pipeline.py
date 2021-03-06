#!/usr/bin/env python

import os
import argparse
import configparser
import csv
import datetime
import drmaa
import glob
import json
import logging
import random
import structlog
import subprocess
import sys
import multiprocessing
import uuid
from pkg_resources import resource_filename
from pprint import pprint

import cpo_pipeline
from cpo_pipeline.logging import now


logger = structlog.get_logger()

def collect_final_outputs(outdir, sample_id):
    final_outputs = {}
    final_outputs['sample_id'] = sample_id
    total_bp_path = os.path.join(
        outdir,
        sample_id,
        'pre-assembly_qc',
        'totalbp'
    )

    try:
        total_bp = cpo_pipeline.assembly.parsers.result_parsers.parse_total_bp(
            total_bp_path
        )
        logger.info(
            "parsed_result_file",
            timestamp=str(now()),
            filename=os.path.abspath(total_bp_path),
        )
    except FileNotFoundError:
        logger.error(
            "output_parsing_failed",
            timestamp=str(now()),
        )
        total_bp = None

    estimated_genome_coverage_stats_path = os.path.join(
        outdir,
        sample_id,
        'pre-assembly_qc',
        'estimated_coverage_stats.tsv'
    )

    try:
        estimated_coverage_stats = cpo_pipeline.assembly.parsers.result_parsers.parse_estimated_coverage_stats(
            estimated_genome_coverage_stats_path
        )
        logger.info(
            "parsed_result_file",
            timestamp=str(now()),
            filename=os.path.abspath(estimated_genome_coverage_stats_path),
        )
    except FileNotFoundError:
        logger.error(
            "output_parsing_failed",
            timestamp=str(now()),
            filename=os.path.abspath(estimated_genome_coverage_stats_path),
        )
        estimated_coverage_stats = {
            'sample_id': sample_id,
            'total_bp': '-',
            'estimated_genome_size': '-',
            'estimated_depth_of_coverage': '-',
        }

    reference_genome_assembly_stats_glob = os.path.join(
        outdir,
        sample_id,
        'reference',
        "*_assembly_stats.txt"
    )
    try:
        [reference_genome_assembly_stats_path] = glob.glob(
            reference_genome_assembly_stats_glob
        )
    except ValueError:
        logger.error(
            "result_parsing_failed",
            timestamp=str(now()),
            filename=str(reference_genome_assembly_stats_glob),
        )
    try:
        reference_genome_assembly_stats = cpo_pipeline.assembly.parsers.result_parsers.parse_reference_genome_assembly_stats(
            reference_genome_assembly_stats_path
        )
    except FileNotFoundError:
        logger.error(
            "output_parsing_failed",
            timestamp=str(now()),
            filename=os.path.abspath(reference_genome_assembly_stats_path),
        )
        reference_genome_assembly_stats = {
            'organism_name': 'Unknown (parsing failed)',
            'infraspecific_name': 'Unknown (parsing failed)',
            'refseq_assembly_accession': 'Unknown (parsing failed)',
            'taxid': 'Unknown (parsing failed)',
            'total_length': 0,
            'contig_count': 0,
            'contig_N50': 0,
        }

    mlst_result_path = os.path.join(
        outdir,
        sample_id,
        'typing',
        'mlst',
        'mlst.tsv'
    )
    try:
        [mlst_result] = cpo_pipeline.typing.parsers.result_parsers.parse_mlst_result(
            mlst_result_path
        )
    except ValueError:
        logger.error(
            "output_parsing_failed",
            timestamp=str(now()),
        )
        mlst_result = {
            'contig_file': os.path.join(
                outdir,
                sample_id,
                'assembly',
                'contigs.fa'
                ),
            'scheme_id': '-',
            'sequence_type': '-',
            'multi_locus_alleles': {
	        'adk': '-',
	        'fumc': '-',
	        'gyrB': '-',
	        'icd': '-',
	        'mdh': '-',
	        'purA': '-',
	        'recA': '-'
            }
        }

    final_outputs['bp'] = total_bp
    final_outputs['est_genome_size'] = estimated_coverage_stats['estimated_genome_size']
    final_outputs['Coverage'] = round(estimated_coverage_stats['estimated_depth_of_coverage'], 2)
    final_outputs['MASH_BEST_HIT'] = " ".join([
        reference_genome_assembly_stats['organism_name'],
        reference_genome_assembly_stats['infraspecific_name'],   
    ])
    final_outputs['MLST_SCHEME'] = mlst_result['scheme_id']
    final_outputs['MLST'] = mlst_result['sequence_type']
    allele_number = 1
    for key, value in mlst_result['multi_locus_alleles'].items():
        final_outputs['MLST_ALLELE_' + str(allele_number)] = key + "(" + value + ")"
        allele_number += 1

    return final_outputs


def main(args):
    """
    """


    config = configparser.ConfigParser()
    config.read(args.config_file)

    analysis_id = uuid.uuid4()

    
    logger.new(
        analysis_id=str(uuid.uuid4()),
        sample_id=args.sample_id,
        pipeline_version=cpo_pipeline.__version__,
    )

    logger.info(
        "analysis_started",
        timestamp=str(now()),
    )

    cpo_pipeline.plasmids.pipeline.main(args)

    cpo_pipeline.assembly.pipeline.main(args)

    cpo_pipeline.typing.pipeline.main(args)

    cpo_pipeline.resistance.pipeline.main(args)

    final_outputs = collect_final_outputs(args.outdir, args.sample_id)
    logger.info(
        "collected_final_outputs",
        final_outputs=final_outputs,
    )

    final_output_path = "/".join([
        args.outdir,
        args.sample_id,
        'final_output.tsv'
    ])

    final_outputs_headers = [
        'sample_id',
        'bp',
        'est_genome_size',
        'Coverage',
        'MASH_BEST_HIT',
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

    with open(final_output_path, 'w+') as f:
        writer = csv.DictWriter(f, fieldnames=final_outputs_headers, delimiter='\t')
        writer.writeheader()
        writer.writerow(final_outputs)


    logger.info(
        "analysis_completed",
        timestamp=str(now()),
    )

def multi(args):
    """
    """
    
    config = configparser.ConfigParser()
    config.read(args.config_file)

    pool = multiprocessing.Pool(int(args.parallel))

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
            parser.add_argument("-e", "--expected-organism-taxid", dest="expected_organism_ncbi_taxid",
                                help="Expected organism NCBI Taxonomy ID")
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

    pool.map(main, arguments)
    
if __name__ == '__main__':
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

    main(args)

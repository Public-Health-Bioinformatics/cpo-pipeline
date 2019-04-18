#!/usr/bin/env python

import argparse
import configparser
import datetime
import logging
import os
import structlog
import sys
import time

from pkg_resources import resource_filename


from cpo_pipeline.drmaa import prepare_job, run_jobs
from cpo_pipeline.logging import now
from cpo_pipeline.typing.parsers import result_parsers
from cpo_pipeline.typing.parsers import input_parsers


logger = structlog.get_logger()

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

    sample_id = args.sample_id
    output_dir = args.outdir
    
    try:
        assembly = args.assembly
    except AttributeError:
        assembly = os.path.join(
            output_dir,
            sample_id,
            'assembly',
            'contigs.fa'
        )

    try:
        mlst_scheme_map_file = args.mlst_scheme_map_file
    except AttributeError:
        mlst_scheme_map_file = resource_filename('data', 'scheme_species_map.tab')
    if not mlst_scheme_map_file:
        mlst_scheme_map_file = resource_filename('data', 'scheme_species_map.tab')

    
    paths = {
        "output_dir": output_dir,
        'logs': os.path.join(
            output_dir,
            sample_id,
            'logs',
        ),
        'mlst_path': os.path.join(
            output_dir,
            sample_id,
            'typing',
            'mlst',
            'mlst.tsv'
        ),
        'mob_recon_path': os.path.join(
            output_dir,
            sample_id,
            'typing',
            'mob_recon'
        ),
        'abricate_plasmidfinder_path': os.path.join(
            output_dir,
            sample_id,
            'typing',
            'abricate',
            'abricate_plasmidfinder.tsv'
        ),
    }

    job_script_path = resource_filename('data', 'job_scripts')

    typing_jobs = [
        {
            'job_name': "_".join(['mlst', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mlst.sh'),
            'args': [
                "--input", assembly,
                "--label", sample_id,
                "--output_file", paths['mlst_path']
            ]
        },
        {
            'job_name': "_".join(['abricate', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'abricate.sh'),
            'args': [
                "--input", assembly,
                "--database", "plasmidfinder",
                "--output_file", paths['abricate_plasmidfinder_path']
            ]
        },
        {
            'job_name': "_".join(['mob_recon', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mob_recon.sh'),
            'args': [
                "--input", assembly,
                "--output_dir", paths['mob_recon_path']
            ]
        }
    ]

    run_jobs(typing_jobs)

    mlst_report = os.path.join(
        output_dir,
        sample_id,
        "typing",
        "mlst",
        "mlst.tsv"
    )
    mlst_hits = result_parsers.parse_mlst_result(mlst_report)
    # TODO: Check that there is only one MLST result in the report, and handle
    #       cases where the report is malformed.
    [mlst_hit] = mlst_hits
    logger.info(
        "parsed_result_file",
        timestamp=str(now()),
        filename=os.path.abspath(mlst_report),
        scheme_id=mlst_hit["scheme_id"],
        sequence_type=mlst_hit["sequence_type"],
    )
    mlst_scheme_map = input_parsers.parse_scheme_species_map(mlst_scheme_map_file)
    mlst_species = "Undefined"
    for scheme in mlst_scheme_map:
        if 'species' in scheme and scheme['scheme_id'] == mlst_hit['scheme_id']:
            mlst_species = scheme['species']

    mob_recon_contig_report_path = os.path.join(
        output_dir,
        sample_id,
        "typing",
        "mob_recon",
        "contig_report.txt"
    )

    mob_recon_contig_report = result_parsers.parse_mob_recon_contig_report(
        mob_recon_contig_report_path
    )
    logger.info(
        "parsed_result_file",
        timestamp=str(now()),
        filename=os.path.abspath(mob_recon_contig_report_path),
        num_records=len(mob_recon_contig_report),
    )
    
    mob_recon_aggregate_report_path = os.path.join(
        output_dir,
        sample_id,
        "typing",
        "mob_recon",
        "mobtyper_aggregate_report.txt"
    )

    mob_recon_aggregate_report = result_parsers.parse_mob_recon_mobtyper_aggregate_report(
        mob_recon_aggregate_report_path
    )
    logger.info(
        "parsed_result_file",
        timestamp=str(now()),
        filename=os.path.abspath(mob_recon_aggregate_report_path),
        num_records=len(mob_recon_aggregate_report),
    )

    def extract_contig_num(contig_id):
        """
        Given a contig_id from a mob_recon contig_report.txt file, return only the contig number.
        Args:
            contig_id (str): contig_id field from mob_recon contig_report.txt
            For example: "contigs.fa|contig00054_len=2672_cov=424.9_corr=0_origname=NODE_54_length_2672_cov_424.949312_pilon_sw=shovill-spades/1.0.1_date=20181024"
        Returns:
            str: contig number.
            For example: "00054"
        """
        prefix = '|contig'
        suffix = '_len='
        prefix_index = contig_id.find(prefix) + len(prefix)
        suffix_index = contig_id.find(suffix)
        contig_num = contig_id[prefix_index:suffix_index]
        return contig_num

    def get_plasmid_contigs(mob_recon_contig_report):
        """
        Given a list of dicts generated by parsing a mob_recon contig_report.txt file,
        return a list of plasmid contigs.
        Args:
            mob_recon_contig_report (list of dict):
        Returns:
            list: plasmid contigs
            For example: ['00021', '00022', '00032', ...]
        """
        plasmid_contigs = []
        for contig_report_record in mob_recon_contig_report:
            contig_num = extract_contig_num(contig_report_record['contig_id'])
            if contig_num not in plasmid_contigs and contig_report_record['rep_type']:
                plasmid_contigs.append(contig_num)
        return plasmid_contigs

    def get_likely_plasmid_contigs(mob_recon_contig_report):
        """
        Given a list of dicts generated by parsing a mob_recon contig_report.txt file,
        return a list of likely plasmid contigs.
        Args:
            mob_recon_contig_report (list of dict):
        Returns:
            list: likely plasmid contigs
            For example: ['00054', '00039', '00061', ...]
        """
        likely_plasmid_contigs = []
        for contig_report_record in mob_recon_contig_report:
            contig_num = extract_contig_num(contig_report_record['contig_id'])
            if contig_num not in likely_plasmid_contigs and not contig_report_record['rep_type']:
                likely_plasmid_contigs.append(contig_num)
        return likely_plasmid_contigs

    def get_plasmid_origins(mob_recon_contig_report):
        """
        Given a list of dicts generated by parsing a mob_recon contig_report.txt file,
        return a list of plasmid origins.
        Args:
            mob_recon_contig_report (list of dict):
        Returns:
            list: plasmid origins
            For example: ['rep_cluster_1254', 'IncL/M', 'IncN', ...]
        """
        origins = []
        for contig_report_record in mob_recon_contig_report:
            if contig_report_record['rep_type']:
                if contig_report_record['rep_type'] not in origins:
                    origins.append(contig_report_record['rep_type'])
        return origins

    plasmid_contigs = get_plasmid_contigs(mob_recon_contig_report)
    likely_plasmid_contigs = get_likely_plasmid_contigs(mob_recon_contig_report)
    origins = get_plasmid_origins(mob_recon_contig_report)


if __name__ == "__main__":
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name)
    parser.add_argument("-i", "--ID", dest="sample_id",
                        help="identifier of the isolate")
    parser.add_argument("-o", "--outdir", dest="outdir", default='./',
                        help="absolute path to output folder")
    parser.add_argument("-a", "--assembly", dest="assembly",
                        help="Path to assembly file.")
    parser.add_argument("--mlst-scheme-map", dest="mlst_scheme_map_file",
                        help="absolute file path to mlst scheme")
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

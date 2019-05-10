#!/usr/bin/env python

import argparse
import os
import datetime
import logging
import structlog
import sys
import time
import configparser
import re

from pkg_resources import resource_filename

from cpo_pipeline.drmaa import prepare_job, run_jobs
from cpo_pipeline.logging import now
from cpo_pipeline.resistance.parsers import result_parsers


logger = structlog.get_logger()

def main(args):
    """
    main entrypoint
    Args:
        args(argparse.Namespace): Parsed command-line arguments.
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
        card_path = args.card_json
    except AttributeError:
        try:
            card_path = config['databases']['card_json']
            if not os.path.exists(card_path):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), card_path
                )
            logger.info(
                "configuration_loaded",
                timestamp=str(now()),
                configuration_attribute="databases/card_json",
                configuration_value=card_path,
            )
        except Exception as e:
            logger.error(
                "configuration_failed",
                timestamp=str(now()),
                configuration_attribute="databases/card_json",
                error_message=str(e),
            )

    try:
        abricate_datadir = args.abricate_datadir
    except AttributeError:
        try:
            abricate_datadir = config['databases']['abricate_datadir']
            if not os.path.exists(abricate_datadir):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), abricate_datadir
                )
            logger.info(
                "configuration_loaded",
                timestamp=str(now()),
                configuration_attribute="databases/abricate_datadir",
                configuration_value=abricate_datadir,
            )
        except Exception as e:
            logger.error(
                "configuration_failed",
                timestamp=str(now()),
                configuration_attribute="databases/abricate_datadir",
                error_message=str(e),
            )

    try:
        abricate_cpo_plasmid_db = args.abricate_cpo_plasmid_db
    except AttributeError:
        try:
            abricate_cpo_plasmid_db = config['databases']['abricate_cpo_plasmid_db']
            if not os.path.exists(os.path.join(abricate_datadir, abricate_cpo_plasmid_db)):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), abricate_cpo_plasmid_db
                )
            logger.info(
                "configuration_loaded",
                timestamp=str(now()),
                configuration_attribute="databases/abricate_cpo_plasmid_db",
                configuration_value=abricate_cpo_plasmid_db,
            )
        except Exception as e:
            logger.error(
                "configuration_failed",
                timestamp=str(now()),
                configuration_attribute="databases/abricate_cpo_plasmid_db",
            )

    paths = {
        "output_dir": output_dir,
        'logs': os.path.join(
            output_dir,
            sample_id,
            'logs',
        ),
        'abricate_path': os.path.join(
            output_dir,
            sample_id,
            'resistance',
            'abricate',
            'abricate.tsv'
        ),
        'rgi_path': os.path.join(
            output_dir,
            sample_id,
            'resistance',
            'rgi'
        ),
    }

    job_script_path = resource_filename('data', 'job_scripts')

    resistance_jobs = [
        {
            'job_name': "_".join(['abricate', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'abricate.sh'),
            'args': [
                "--input", assembly,
                "--datadir", abricate_datadir,
                "--database", abricate_cpo_plasmid_db,
                "--output_file", paths['abricate_path']
            ]
        },
        {
            'job_name': "_".join(['rgi', sample_id]),
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'rgi.sh'),
            'args': [
                "--input", assembly,
                "--card_json", card_path,
                "--output_dir", paths['rgi_path']
            ]
        }
    ]

    run_jobs(resistance_jobs)

    abricate_report_path = os.path.join(
        output_dir,
        sample_id,
        "resistance",
        "abricate",
        "abricate.tsv"
    )
    abricate_report = result_parsers.parse_abricate_result(abricate_report_path)
    logger.info(
        "parsed_result_file",
        timestamp=str(datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat()),
        filename=os.path.abspath(abricate_report_path),
        resistance_genes=[
            {
                key: record[key] for key in [
                    "gene",
                    "accession",
                    "database",
                    "percent_coverage",
                    "percent_identity",
                ]
            }
            for record in abricate_report
        ]
    )
    
    rgi_report_path = os.path.join(
        output_dir,
        sample_id,
        "resistance",
        "rgi",
        "rgi.txt"
    )
    rgi_report = result_parsers.parse_rgi_result_txt(rgi_report_path)
    logger.info(
        "parsed_result_file",
        timestamp=str(datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat()),
        filename=os.path.abspath(rgi_report_path),
        resistance_genes=[
            {
                key: record[key] for key in [
                    "best_hit_aro",
                    "aro",
                ]
            }
            for record in rgi_report
        ]
    )
    
    def get_abricate_carbapenemases(abricate_report):
        """
        Given a list of dicts generated by parsing an abricate report file,
        return a list of carbapenemases.
        Args:
            abricate_report (list of dict):
        Returns:
            list: likely plasmid contigs
            For example: ['NDM-1', '', '', ...]
        """
        abricate_carbapenemases = []
        for abricate_report_record in abricate_report:
            abricate_carbapenemases.append(abricate_report_record['gene'])
        return abricate_carbapenemases

    def get_rgi_carbapenemases(rgi_report):
        """
        Given a list of dicts generated by parsing an rgi report file,
        return a list of carbapenemases.
        Args:
            rgi_report (list of dict):
        Returns:
            list: likely plasmid contigs
            For example: ['', '', '', ...]
        """
        rgi_carbapenemases = []
        for rgi_report_record in rgi_report:
            if re.search("carbapenem", rgi_report_record['drug_class']):
                rgi_carbapenemases.append(rgi_report_record['best_hit_aro'])
        return rgi_carbapenemases

if __name__ == "__main__":
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name)
    parser.add_argument("-i", "--ID", dest="sample_id",
                        help="identifier of the isolate", required=True)
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="absolute path to output folder", required=True)
    parser.add_argument("-a", "--assembly", dest="assembly",
                        help="Path to assembly file.", required=True)
    parser.add_argument("--card-json", dest="card_json",
                        help="absolute path to card database (json format)")
    parser.add_argument("--abricate-cpo-plamid-db", dest="abricate_cpo_plasmid_db",
                        help="absolute path to card database (json format)")
    parser.add_argument("--abricate-datadir", dest="abricate_datadir",
                        help="absolute path to card database (json format)")
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

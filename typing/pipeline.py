#!/usr/bin/env python

'''
This script is a wrapper for module two, part 1: typing from assemblies. 

It uses mlst, mobsuite for sequence typing. 

Example usage:

  pipeline.py --id BC11-Kpn005 --assembly BC11-Kpn005_S2.fa --output output --expected-species "Klebsiella"

Requires the pipeline_prediction.sh script. the directory of where it's store can be specified using -k.
'''

import subprocess
import optparse
import os
import datetime
import sys
import time
import urllib.request
import gzip
import collections
import configparser
import drmaa

from parsers import result_parsers
from parsers import input_parsers


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

def main():

    config = configparser.ConfigParser()
    config.read(os.path.dirname(os.path.realpath(sys.argv[0])) + '/config.ini')

    parser = optparse.OptionParser("Usage: %prog [options] arg1 arg2 ...")
    #required
    parser.add_option("-i", "--id", dest="id", type="string", help="identifier of the isolate")
    parser.add_option("-a", "--assembly", dest="assembly", type="string", help="Path to assembly file.")
    parser.add_option("-o", "--output", dest="output", default='./', type="string", help="absolute path to output folder")    

    parser.add_option("-s", "--mlst-scheme-map", dest="mlst_scheme_map_file", default=config['databases']['mlst-scheme-map'], type="string", help="absolute file path to mlst scheme")
    parser.add_option("-k", "--script-path", dest="script_path", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")
    
    
    (options, args) = parser.parse_args()
    curDir = os.getcwd()
    ID = str(options.id).lstrip().rstrip()
    assembly = options.assembly
    script_path = options.script_path
    mlst_scheme_map_file = options.mlst_scheme_map_file
    outputDir = options.output


    print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nAssembly: " + assembly)

    file_paths = {
        'mlst_path': '/'.join([outputDir, ID, 'typing', 'mlst', 'mlst.tsv']),
        'mob_recon_path': '/'.join([outputDir, ID, 'typing', 'mob_recon']),
    }
    
    job_script_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),  'job_scripts')
    
    typing_jobs = [
        {
            'job_name': 'mlst',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mlst.sh'),
            'args': [
                "--input", assembly,
                "--output_file", file_paths['mlst_path']
            ]
        },
        {
            'job_name': 'mob_recon',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mob_recon.sh'),
            'args': [
                "--input", assembly,
                "--output_dir", file_paths['mob_recon_path']
            ]
        }
    ]

    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in typing_jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)    
    
    print("identifying MLST")
    mlst_report = "/".join([outputDir, ID, "typing", "mlst", "mlst.tsv"]) 
    mlstHits = result_parsers.parse_mlst_result(mlst_report)
    # TODO: Check that there is only one MLST result in the report, and handle
    #       cases where the report is malformed.
    mlstHit = mlstHits[0]
    mlst_scheme_map = input_parsers.parse_scheme_species_map(mlst_scheme_map_file)
    mlst_species = "Undefined"
    for scheme in mlst_scheme_map:
        if 'species' in scheme and scheme['scheme_id'] == mlstHit['scheme_id']:
            mlst_species = scheme['species']

    print("identifying plasmid contigs and amr genes")

    mob_recon_contig_report_path = "/".join([outputDir, ID, "typing", "mob_recon", "contig_report.txt"])
    mob_recon_contig_report = result_parsers.parse_mob_recon_contig_report(mob_recon_contig_report_path)

    mob_recon_aggregate_report_path = "/".join([outputDir, ID, "typing", "mob_recon", "mobtyper_aggregate_report.txt"])
    mob_recon_aggregate_report = result_parsers.parse_mob_recon_mobtyper_aggregate_report(mob_recon_aggregate_report_path)
    

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

    def get_likely_plasmid_contigs(mob_recon_report):
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
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

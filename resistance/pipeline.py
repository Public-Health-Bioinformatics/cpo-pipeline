#!/usr/bin/env python

'''
This script is a wrapper for resistance gene idenfification from assemblies. 

It uses abricate and rgi for AMR profile prediction and plasmid predictions. 

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
    parser.add_option("-k", "--script-path", dest="script_path", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")
    parser.add_option("-c", "--card-path", dest="card_path", default=config['databases']['card'], type="string", help="absolute file path to card.json db")
    parser.add_option("-d", "--abricate-datadir", dest="abricate_datadir", default=config['databases']['abricate-datadir'], type="string", help="absolute file path to directory where abricate dbs are stored")
    parser.add_option("-p", "--abricate-cpo-plasmid-db", dest="abricate_cpo_plasmid_db", default=config['databases']['abricate-cpo-plasmid-db'], type="string", help="name of abricate cpo plasmid db to use")
    
    
    (options, args) = parser.parse_args()
    curDir = os.getcwd()
    ID = str(options.id).lstrip().rstrip()
    assembly = options.assembly
    script_path = options.script_path
    card_path = options.card_path
    abricate_datadir = options.abricate_datadir
    abricate_cpo_plasmid_db = options.abricate_cpo_plasmid_db
    outputDir = options.output


    print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nAssembly: " + assembly)

    file_paths = {
        'abricate_path': '/'.join([outputDir, ID, 'resistance', 'abricate', 'abricate.tsv']),
        'rgi_path': "/".join([outputDir, ID, 'resistance', 'rgi', 'rgi'])
    }
    
    job_script_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),  'job_scripts')
    
    resistance_jobs = [
        {
            'job_name': 'abricate',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'abricate.sh'),
            'args': [
                "--input", assembly,
                "--datadir", abricate_datadir,
                "--database", abricate_cpo_plasmid_db,
                "--output_file", file_paths['abricate_path']
            ]
        },
        {
            'job_name': 'rgi',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'rgi.sh'),
            'args': [
                "--input", assembly,
                "--card_json", card_path,
                "--output_file", file_paths['rgi_path']
            ]
        }
    ]
    
    with drmaa.Session() as session:
        prepared_jobs = [prepare_job(job, session) for job in resistance_jobs]
        running_jobs = [session.runJob(job) for job in prepared_jobs]
        for job_id in running_jobs:
            print('Your job has been submitted with ID %s' % job_id)
        session.synchronize(running_jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    
    abricate_report_path = "/".join([outputDir, ID, "resistance", "abricate", "abricate.tsv"])
    abricate_report = result_parsers.parse_abricate_result(abricate_report_path)
    
    rgi_report_path = "/".join([outputDir, ID, "resistance", "rgi", "rgi.txt"])
    rgi_report = result_parsers.parse_rgi_result_txt(rgi_report_path)

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
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

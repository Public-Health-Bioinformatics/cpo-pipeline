#!/usr/bin/env python

'''
This script is a wrapper for module two, part 1: typing and resistance gene idenfification from assemblies. 

It uses mlst, resfinder, plasmidfinder, rgi, mobsuite for sequence typing, AMR profile prediction and plasmid predictions. 
The output is a TSV summary of the results. 

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

from parsers import result_parsers
from parsers import input_parsers


def execute(command, curDir):
    process = subprocess.Popen(command, shell=False, cwd=curDir, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        raise subprocess.CalledProcessError(exitCode, command)


def main():

    config = configparser.ConfigParser()
    config.read(os.path.dirname(os.path.realpath(sys.argv[0])) + '/config.ini')

    #parses some parameters
    parser = optparse.OptionParser("Usage: %prog [options] arg1 arg2 ...")
    #required
    #MLSTHIT, mobsuite, resfinder, rgi, mlstscheme
    parser.add_option("-i", "--id", dest="id", type="string", help="identifier of the isolate")
    parser.add_option("-a", "--assembly", dest="assembly", type="string", help="Path to assembly file.")
    parser.add_option("-o", "--output", dest="output", default='./', type="string", help="absolute path to output folder")    

    parser.add_option("-s", "--mlst-scheme-map", dest="mlst_scheme_map_file", default=config['databases']['mlst-scheme-map'], type="string", help="absolute file path to mlst scheme")
    parser.add_option("-k", "--script-path", dest="script_path", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")
    parser.add_option("-c", "--card-path", dest="card_path", default=config['databases']['card'], type="string", help="absolute file path to card.json db")
    parser.add_option("-d", "--abricate-datadir", dest="abricate_datadir", default=config['databases']['abricate-datadir'], type="string", help="absolute file path to directory where abricate dbs are stored")
    parser.add_option("-p", "--abricate-cpo-plasmid-db", dest="abricate_cpo_plasmid_db", default=config['databases']['abricate-cpo-plasmid-db'], type="string", help="name of abricate cpo plasmid db to use")
    parser.add_option("-e", "--expected-species", dest="expected_species", default="NA/NA/NA", type="string", help="expected species of the isolate")
    
    # parallelization, useless, these are hard coded to 8cores/64G RAM
    # parser.add_option("-t", "--threads", dest="threads", default=8, type="int", help="number of cpu to use")
    # parser.add_option("-p", "--memory", dest="memory", default=64, type="int", help="memory to use in GB")
    
    (options, args) = parser.parse_args()
    curDir = os.getcwd()
    ID = str(options.id).lstrip().rstrip()
    assembly = options.assembly
    expected_species = options.expected_species
    script_path = options.script_path
    cardPath = options.card_path
    abricate_datadir = options.abricate_datadir
    abricate_cpo_plasmid_db = options.abricate_cpo_plasmid_db
    mlst_scheme_map_file = options.mlst_scheme_map_file
    outputDir = options.output


    print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nAssembly: " + assembly)

    print("running pipeline_typing.sh")
    #input parameters: 1=id, 2=assembly, 3=output, 4=cardPath, 5=abricate_datadir, 6=abricate_cpo_plasmid_db
    cmd = [script_path + "/pipeline_typing.sh", ID, assembly, outputDir, cardPath, abricate_datadir, abricate_cpo_plasmid_db]
    result = execute(cmd, curDir)
    
    print("step 3: parsing mlst, plasmid, and amr results")
    
    print("identifying MLST")
    mlst_report = outputDir + "/typing/" + ID + "/" + ID + ".mlst/" + ID + ".mlst" 
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

    mob_recon_contig_report_path = outputDir + "/typing/" + ID + "/" + ID + ".recon/" + "contig_report.txt" 
    mob_recon_contig_report = result_parsers.parse_mob_recon_contig_report(mob_recon_contig_report_path)

    mob_recon_aggregate_report_path = outputDir + "/typing/" + ID + "/" + ID + ".recon/" + "mobtyper_aggregate_report.txt"
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
    
    abricate_report_path = outputDir + "/resistance/" + ID + "/" + ID + ".cp"
    abricate_report = result_parsers.parse_abricate_result(abricate_report_path)
    
    rgi_report_path = outputDir + "/resistance/" + ID + "/" + ID + ".rgi.txt"
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

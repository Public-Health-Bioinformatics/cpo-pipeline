#!/usr/bin/env python

'''
This script is a wrapper for module one of the cpo-pipeline including QC and Assembly.
It uses Mash2.0 and fastqc to check for sequence contamination, quality information and identify a reference genome.
Then attempts to assemble the reads, attempting to filter contamination away if required.

Example usage:

  pipeline.py -i BC11-Kpn005_S2 --R1 BC11-Kpn005_S2_L001_R1_001.fastq.gz --R2 BC11-Kpn005_S2_L001_R2_001.fastq.gz -o output_dir

Requires pipeline_qc.sh, pipeline_assembly.sh, pipeline_assembly_contaminant.sh. where these scripts are located can be specified with -k. 
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
import json
import configparser
import pprint
import shutil
import errno
import re
import drmaa
import glob

from parsers import result_parsers

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
    if (complete_single / total) >= busco_complete_single_cutoff and \
       (complete_duplicate / total) <= busco_complete_duplicate_cutoff:
        return True
    else:
        return False
    
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
    if total_length <= (reference_length * (1 + assembly_length_cutoff)) and \
       total_length >= (reference_length * (1 - assembly_length_cutoff)) and \
       assembly_percent_gc <= (reference_percent_gc * (1 + percent_gc_cutoff)) and \
       assembly_percent_gc >= (reference_percent_gc * (1 - percent_gc_cutoff)) and \
       genome_fraction_percent >= genome_fraction_percent_cutoff:
        return True
    else:
        return False

def fastqc_qc_check(fastqc_results):
    """
    Args:
        fastqc_results (dict): FastQC results
    Returns:
        boolean: Sequence data passes our FastQC quality criteria
    """
    if fastqc_results["basic_statistics"] == "PASS" and \
       fastqc_results["per_base_sequence_quality"] == "PASS" and \
       fastqc_results["sequence_length_distribution"] == "PASS":
        return True
    else:
        return False

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
        urllib.request.urlretrieve(assembly_stat_url, "/".join([download_path, assembly + "_assembly_stats.txt"]))

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
    
    #parses some parameters
    parser = optparse.OptionParser("Usage: %prog [options] arg1 arg2 ...")
    parser.add_option("-i", "--id", dest="id", type="string", help="identifier of the isolate")    
    parser.add_option("-1", "--R1", dest="R1", type="string", help="absolute file path forward read (R1)")
    parser.add_option("-2", "--R2", dest="R2", type="string", help="absolute file path to reverse read (R2)")
    parser.add_option("-g", "--mash-genomedb", dest="mashGenomeRefDB", default = config['databases']['mash-genomedb'], type="string", help="absolute path to mash reference database")
    parser.add_option("-p", "--mash-plasmiddb", dest="mashPlasmidRefDB", default = config['databases']['mash-plasmiddb'], type="string", help="absolute path to mash reference database")
    parser.add_option("-b", "--busco-db", dest="buscodb", default = config['databases']['busco-db'], type="string", help="absolute path to busco reference database")

    parser.add_option("-o", "--output", dest="output", default='./', type="string", help="absolute path to output folder")
    parser.add_option("-s", "--script-path", dest="script_path", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")
    

    (options,args) = parser.parse_args()

    curDir = os.getcwd()
    outputDir = options.output
    mash_genome_db = options.mashGenomeRefDB
    mash_plasmid_db=options.mashPlasmidRefDB
    script_path = options.script_path
    buscodb = options.buscodb
    ID = options.id
    R1 = options.R1
    R2 = options.R2

    
    notes = []
    #init the output list
    output = []
    
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
        "mash_hits_genome_score_cutoff":300, #genome mash will include all hits with scores (top hit score - $thisvalue)
        "mash_hits_plasmid_score_cutoff":100, #plasmid mash will include all hits with scores (top hit score - $thisvalue)
        "coverage_cutoff":30, #sequencing coverage greater than ($thisvalue) will pass the QC
        "quast_assembly_length_cutoff":0.10, #QUAST QC: assembly length within +-($thisvalue) percent in reference to reference length will pass the QC 
        "quast_percent_gc_cutoff":0.05, #QUAST QC: percent GC within +-($thisvalue) percent in reference to reference percent GC will pass the QC 
        "quast_genome_fraction_percent_cutoff":0.90, #QUAST QC: genome_fraction_percent greater than ($thisvalue) will pass the QC
        "busco_complete_single_cutoff":0.90, #BUSCO QC: complete single genes greater than ($thisvalue) percent will pass the QC
        "busco_complete_duplicate_cutoff":0.10 #BUSCO QC: complete duplicate genes less than ($thisvalue) percent will pass the QC
    }
    
    file_paths = {
        "mash_genome_path": "/".join([outputDir, ID, "pre-assembly_qc", "mashscreen.genome.tsv"]),
        "mash_plasmid_path": "/".join([outputDir, ID, "pre-assembly_qc", "mashscreen.plasmid.tsv"]),
        "fastqc_output_path": "/".join([outputDir, ID, "pre-assembly_qc", "fastqc"]),
        "totalbp_path": "/".join([outputDir, ID, "pre-assembly_qc", "totalbp"]),
        "reference_genome_path": "/".join([outputDir, ID, "reference"]),
        "assembly_path": "/".join([outputDir, ID, "assembly"]),
        "busco_path": "/".join([outputDir, ID, "post-assembly_qc", "busco"]),
        "quast_path": "/".join([outputDir, ID, "post-assembly_qc", "quast"]),
    }
    
    
    print(str(datetime.datetime.now()))
    print("ID: " + ID)
    print("R1: " + R1)
    print("R2: " + R2)
    
    print("step 1: preassembly QC")

    job_script_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),  'job_scripts')

    pre_assembly_qc_jobs = [
        {
            'job_name': 'mash_screen',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_screen.sh'),
            'args': [
                "--R1", R1,
                "--R2", R2,
                "--queries", mash_genome_db,
                "--output_file", file_paths['mash_genome_path']
            ],
        },
        {
            'job_name': 'mash_screen',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'mash_screen.sh'),
            'args': [
                "--R1", R1,
                "--R2", R2,
                "--queries", mash_plasmid_db,
                "--output_file", file_paths['mash_plasmid_path']
            ],
        },
        {
            'job_name': 'fastqc',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'fastqc.sh'),
            'args': [
                "--R1", R1,
                "--R2", R2,
                "--output_dir", file_paths['fastqc_output_path']
            ],
        },
        {
            'job_name': 'seqtk_totalbp',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'seqtk_totalbp.sh'),
            'args': [
                "--R1", R1,
                "--R2", R2,
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
    mash_hits_score_threshold = int(mash_hits[0]['shared_hashes'].split("/")[0]) - int(qc_thresholds["mash_hits_genome_score_cutoff"])

    def score_above_threshold(mash_result, score_threshold):
        score = int(mash_result['shared_hashes'].split("/")[0])
        if (score >= score_threshold and mash_result['query_comment'].find("phiX") == -1):
            return True
        else:
            return False
        
    filtered_mash_hits = list(filter(
        lambda x: score_above_threshold(x, mash_hits_score_threshold),
        mash_hits))
    
    # parse plasmid mash
    mash_plasmid_hits = result_parsers.parse_mash_result(file_paths["mash_plasmid_path"])
    mash_plasmid_hits = sorted(mash_plasmid_hits, key=lambda k: k['identity'], reverse=True)
    # 'shared_hashes' field is string in format '935/1000'
    # Threshold is 100 below highest numerator (ie. '935/100' -> 835)
    mash_plasmid_hits_score_threshold = int(mash_plasmid_hits[0]['shared_hashes'].split("/")[0]) - int(qc_thresholds["mash_hits_plasmid_score_cutoff"])
    filtered_mash_plasmid_hits = list(filter(
        lambda x: score_above_threshold(x, mash_plasmid_hits_score_threshold),
        mash_plasmid_hits))
    
    # parse fastqc
    fastqc_R1 = result_parsers.parse_fastqc_result(file_paths['fastqc_output_path'] + "/" + ID + "_R1_fastqc/summary.txt") 
    fastqc_R2 = result_parsers.parse_fastqc_result(file_paths['fastqc_output_path'] + "/" + ID + "_R2_fastqc/summary.txt") 

    #all the qC result are parsed now, lets do some QC logic
    #look at mash results first
    if (len(filtered_mash_hits) > 1):
        qc_verdicts["multiple_species_contamination"] = True
    else:
        qc_verdicts["multiple_species_contamination"] = False
    
    
    if (len(filtered_mash_plasmid_hits) > 0):
        qc_verdicts["fastq_contains_plasmids"] = True
    else:
        qc_verdicts["fastq_contains_plasmids"] = False

    #look at fastqc results
    qc_verdicts["acceptable_fastqc_forward"] = fastqc_qc_check(fastqc_R1) 
    qc_verdicts["acceptable_fastqc_reverse"] = fastqc_qc_check(fastqc_R2) 
    
    #download a reference genome
    print("Downloading reference genomes")
        
    reference_genomes = []
    # build the save paths
    try:
        os.makedirs(file_paths['reference_genome_path'])
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    if (not qc_verdicts["multiple_species_contamination"]):
        for mash_hit in filtered_mash_hits:
            if not re.match("phiX", mash_hit['query_comment']):
                download_mash_hit(mash_hit, file_paths['reference_genome_path'])
                reference_genomes.append(mash_hit['query_id'])
        
    else: #throw an error if it contains contaminations
        print("Contaminated Genome assembly...resequencing required")
        raise Exception("contamination and mislabeling...crashing")
        
    #check to make sure we ONLY have ONE reference.
    if (len(reference_genomes) > 1 ):
        raise Exception ("there are multiple reference genomes")
    elif (len(reference_genomes) == 0):
        raise Exception ("no reference genome identified")
    
    # now we estimate our coverage using total reads and expected genome size
    expected_genome_size = result_parsers.parse_reference_genome_stats(glob.glob(file_paths["reference_genome_path"] + "/*_assembly_stats.txt")[0])
    total_bp = result_parsers.parse_total_bp(file_paths["totalbp_path"])
    
    coverage = total_bp / expected_genome_size
            
    if (coverage >= int(qc_thresholds["coverage_cutoff"])):
        qc_verdicts["acceptable_coverage"] = True

    
    print("step 2: genome assembly and QC")

    assembly_jobs = [
        {
            'job_name': 'shovill',
            'native_specification': '-pe smp 16 -l h_vmem=4G',
            'remote_command': os.path.join(job_script_path, 'shovill.sh'),
            'args': [
                "--R1", R1,
                "--R2", R2,
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
                "--database", buscodb,
                "--output_dir", file_paths['busco_path']
            ]
        },
        {
            'job_name': 'quast',
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'quast.sh'),
            'args': [
                "--input", "/".join([file_paths['assembly_path'], "contigs.fa"]),
                "--reference_genome", "/".join([file_paths['reference_genome_path'], reference_genomes[0]]),
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
    busco_results = result_parsers.parse_busco_result(file_paths["busco_path"] + "/run_busco/short_summary_busco.txt")
    quast_results = result_parsers.parse_quast_result(file_paths["quast_path"] + "/report.txt")

    qc_verdicts["acceptable_busco_assembly_metrics"] = busco_qc_check(busco_results, qc_thresholds)
    qc_verdicts["acceptable_quast_assembly_metrics"] = quast_qc_check(quast_results, qc_thresholds)

    #print QC results to screen
    print("total bases: " + str(total_bp))
    print("expected genome size: " + str(expected_genome_size))
    print("coverage: " + str(coverage))
    print("")
    for key,value in qc_verdicts.items():
        print (str(key) + ": " + str(value))        

    
if __name__ == "__main__":
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

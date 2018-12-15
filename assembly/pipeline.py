#!/usr/bin/env python

'''
This script is a wrapper for module one of the cpo-pipeline including QC and Assembly.
It uses Mash2.0, Kraken2.0 and fastqc to check for sequence contamination, quality information and identify a reference genome.
Then attempts to assemble the reads, attempting to filter contamination away if required.

Example usage:

  pipeline.py -i BC11-Kpn005_S2 -f BC11-Kpn005_S2_L001_R1_001.fastq.gz -r BC11-Kpn005_S2_L001_R2_001.fastq.gz -o output_dir -e "Klebsiella pneumoniae" 

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

from parsers import result_parsers

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


def busco_results_qc_check(busco_results):
    '''
    BUSCO PASS CRITERIA:
    1. complete singles > 90% of total genes
    2. complte duplicates < 90% of total genes
    Args:
        busco_results (dict): Busco results
    Returns:
        boolean: Assembly passes our BUSCO quality criteria
    '''
    complete_single = busco_results['complete_single']
    total = busco_results['total']
    busco_complete_single_cutoff = busco_results['busco_complete_single_cutoff']
    busco_complete_duplicate_cutoff = busco_results['busco_complete_duplicate_cutoff']
    if (complete_single / total) >= busco_complete_single_cutoff and \
       (complete_duplicate / total) <= busco_complete_duplicate_cutoff:
        return True
    else:
        return False
    
def quast_results_qc_check(quast_results):
    '''
    QUAST PASS CRITERIA:
    1. total length vs reference length +-10%
    2. percent gc versus reference percent gc +- 5%
    3. genome fraction percent > 90
    Args:
        quast_results (dict): Quast results

    Returns:
        boolean: Assembly passes our QUAST quality criteria
    '''
    total_length = quast_results['total_length']
    reference_length = quast_results['reference_length']
    assembly_length_cutoff = quast_results['quast_assembly_length_cutoff']
    assembly_percent_gc = quast_results['percent_GC']
    reference_percent_gc = quast_results['reference_percent_GC']
    percent_gc_cutoff = quast_results['quast_percent_gc_cutoff']
    genome_fraction_percent = quast_results['genome_fraction_percent']
    genome_fraction_percent_cutoff = quast_results['genome_fraction_percent_cutoff']
    if total_length <= (reference_length * (1 + assembly_length_cutoff)) and \
       total_length >= (reference_length * (1 - assembly_length_cutoff)) and \
       assembly_percent_gc <= (reference_percent_gc * (1 + percent_gc_cutoff)) and \
       assembly_percent_gc >= (reference_percent_gc * (1 - quast_percent_gc_cutoff)) and \
       genome_fraction_percent >= genome_fraction_percent_cutoff:
        return True
    else:
        return False

def mash_query_id_to_ncbi_ftp_path(query_id):
    """
    Args:
        query_id (str): Mash query ID (column 5 of mash screen report)
    Returns:
        list: Directory names used to locate reference genome on ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
        For example:
            "GCF/001/022/155"
    """
    prefix = query_id.split('_')[0]
    digits = query_id.split('_')[1].split('.')[0]
    path_list = [prefix] + [digits[i:i+3] for i in range(0, len(digits), 3)]
    
    return "/".join(path_list)

    
def main():
    
    config = configparser.ConfigParser()
    config.read(os.path.dirname(os.path.realpath(sys.argv[0])) + '/config.ini')
    
    #parses some parameters
    parser = optparse.OptionParser("Usage: %prog [options] arg1 arg2 ...")
    parser.add_option("-i", "--id", dest="id", type="string", help="identifier of the isolate")    
    parser.add_option("-f", "--forward", dest="R1", type="string", help="absolute file path forward read (R1)")
    parser.add_option("-r", "--reverse", dest="R2", type="string", help="absolute file path to reverse read (R2)")
    parser.add_option("-m", "--mash-genomedb", dest="mashGenomeRefDB", default = config['databases']['mash-genomedb'], type="string", help="absolute path to mash reference database")
    parser.add_option("-n", "--mash-plasmiddb", dest="mashPlasmidRefDB", default = config['databases']['mash-plasmiddb'], type="string", help="absolute path to mash reference database")
    parser.add_option("-z", "--kraken2-genomedb", dest="kraken2GenomeRefDB", default = config['databases']['kraken2-genomedb'], type="string", help="absolute path to kraken reference database")
    parser.add_option("-v", "--kraken2-plasmiddb", dest="kraken2PlasmidRefDB", default = config['databases']['kraken2-plasmiddb'], type="string", help="absolute path to kraken reference database")
    parser.add_option("-x", "--busco-db", dest="buscodb", default = config['databases']['busco-db'], type="string", help="absolute path to busco reference database")

    parser.add_option("-o", "--output", dest="output", default='./', type="string", help="absolute path to output folder")
    parser.add_option("-k", "--script-path", dest="script_path", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")

    #used for parsing 
    parser.add_option("-e", "--expected", dest="expectedSpecies", default="NA/NA/NA", type="string", help="expected species of the isolate")
    

    (options,args) = parser.parse_args()

    curDir = os.getcwd()
    outputDir = options.output
    mashdb = options.mashGenomeRefDB
    mashplasmiddb=options.mashPlasmidRefDB
    kraken2db = options.kraken2GenomeRefDB
    kraken2plasmiddb=options.kraken2PlasmidRefDB
    expectedSpecies = options.expectedSpecies
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
        "same_as_expected_species": None,
        "fastq_contains_plasmids": None,
        "acceptable_coverage": None,
        "acceptable_fastqc_forward": None,
        "acceptable_fastqc_reverse": None,
        "acceptable_quast_assembly_metrics": None,
        "acceptable_busco_assembly_metrics": None
    }
    
    qc_cutoffs = {
        "mash_hits_genome_score_cutoff":300, #genome mash will include all hits with scores (top hit score - $thisvalue)
        "mash_hits_plasmid_score_cutoff":100, #plasmid mash will include all hits with scores (top hit score - $thisvalue)
        "coverage_cutoff":30, #sequencing coverage greater than ($thisvalue) will pass the QC
        "quast_assembly_length_cutoff":0.10, #QUAST QC: assembly length within +-($thisvalue) percent in reference to reference length will pass the QC 
        "quast_percent_gc_cutoff":0.05, #QUAST QC: percent GC within +-($thisvalue) percent in reference to reference percent GC will pass the QC 
        "genome_fraction_percent_cutoff":0.90, #QUAST QC: genome_fraction_percent greater than ($thisvalue) will pass the QC
        "busco_complete_single_cutoff":0.90, #BUSCO QC: complete single genes greater than ($thisvalue) percent will pass the QC
        "busco_complete_duplicate_cutoff":0.10 #BUSCO QC: complete duplicate genes less than ($thisvalue) percent will pass the QC
    }
    
    file_paths = {
        "mash_genome_path": "/".join([outputDir, ID, "pre-assembly_qc", "mashscreen.genome.tsv"]),
        "mash_plasmid_path": "/".join([outputDir, ID, "pre-assembly_qc", "mashscreen.plasmid.tsv"]),
        "fastqc_output_path": "/".join([outputDir, ID, "pre-assembly_qc", "fastqc"]),
        "totalbp_path": "/".join([outputDir, ID, "pre-assembly_qc", "totalbp"]),
        "reference_genome_fasta_path": (""), #built later
        "reference_genome_stat_path": (""), #built later
        "busco_path": "/".join([outputDir, ID, "post-assembly_qc", "busco"]),
        "quast_path": "/".join([outputDir, ID, "post-assembly_qc", "quast"]),
    }
    
    
    print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nR1: " + R1 + "\nR2: " + R2)
    output.append(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nR1: " + R1 + "\nR2: " + R2)
    
    print("step 1: preassembly QC")
    
    print("running mash_screen.sh on genomic db")
    cmd = [script_path + "/job_scripts/mash_screen.sh", "--R1", R1, "--R2", R2,
           "--queries", mashdb, "--output_file", file_paths['mash_genome_path']]
    _ = execute(cmd, curDir)
    
    print("running mash_screen.sh on plasmid db")
    cmd = [script_path + "/job_scripts/mash_screen.sh", "--R1", R1, "--R2", R2,
           "--queries", mashplasmiddb, "--output_file", file_paths['mash_plasmid_path']]
    _ = execute(cmd, curDir)
    
    print("running fastqc")
    cmd = [script_path + "/job_scripts/fastqc.sh", "--R1", R1, "--R2", R2,
           "--output_dir", file_paths['fastqc_output_path']]
    _ = execute(cmd, curDir)
    
    print("running seqtk to calculate totalbp")
    cmd = [script_path + "/job_scripts/seqtk_totalbp.sh", "--R1", R1, "--R2", R2,
           "--output_file", file_paths['totalbp_path']]
    _ = execute(cmd, curDir)
    
    print("Parsing the QC results")
    #parse genome mash results
    mash_hits = result_parsers.parse_mash_result(file_paths["mash_genome_path"])
    mash_hits = sorted(mash_hits, key=lambda k: k['identity'], reverse=True)
    # 'shared_hashes' field is string in format '935/1000'
    # Threshold is 300 below highest numerator (ie. '935/100' -> 635)
    mash_hits_score_threshold = int(mash_hits[0]['shared_hashes'].split("/")[0]) - int(qc_cutoffs["mash_hits_genome_score_cutoff"])
    print("*** mash_hits_score_threshold: " + str(mash_hits_score_threshold))
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
    mash_plasmid_hits_score_threshold = int(mash_plasmid_hits[0]['shared_hashes'].split("/")[0]) - int(qc_cutoffs["mash_hits_plasmid_score_cutoff"])
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
    
    for mash_hit in filtered_mash_hits:
        species = mash_hit['query_comment']
        if (species.find(expectedSpecies) > -1):
            qc_verdicts["same_as_expected_species"] = True
    
    if (len(filtered_mash_plasmid_hits) > 0):
        qc_verdicts["fastq_contains_plasmids"] = True

    #look at fastqc results
    if (fastqc_R1["basic_statistics"] == "PASS" and fastqc_R1["per_base_sequence_quality"] == "PASS" and fastqc_R1["sequence_length_distribution"] == "PASS" ):
        qc_verdicts["acceptable_fastqc_forward"] = True 
    if (fastqc_R2["basic_statistics"] == "PASS" and fastqc_R2["per_base_sequence_quality"] == "PASS" and fastqc_R2["sequence_length_distribution"] == "PASS" ):
        qc_verdicts["acceptable_fastqc_reverse"] = True 
    pprint.pprint(qc_verdicts)
    #download a reference genome
    print("Downloading reference genomes")
    reference_genomes = []
    if (not qc_verdicts["multiple_species_contamination"]):
        for mash_hit in filtered_mash_hits: #for all the mash hits, aka reference genomes
            query_id = mash_hit['query_id'] #hit genome within mash results
            species_name_start = int(mash_hit['query_comment'].index(".")) + 3 #find the start of species name within query_comment column
            species_name_stop = int (mash_hit['query_comment'].index(",")) #find the end of the species name within query_comment column
            if (mash_hit['query_comment'].find("phiX") > -1):
                species = "PhiX" #phix
            else:
                species = str(mash_hit['query_comment'])[species_name_start: species_name_stop] #assign proper species name for reference genome file name
                # find gcf accession
                # TODO: document this or clean it up to be more readable

                ncbi_ftp_path = mash_query_id_to_ncbi_ftp_path(query_id)

                assembly = query_id[:query_id.find("_genomic.fna.gz")] #find the assembly name

                #build the urls
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
                
                # build the save paths
                reference_dir_path = os.path.abspath("/".join([outputDir, ID, "reference"]))
                try:
                    os.makedirs(reference_dir_path)
                except OSError as exc:
                    if exc.errno != errno.EEXIST:
                        raise
                    pass
                
                file_paths["reference_genome_fasta_path"] = reference_dir_path + "/" + species.replace(" ","_") + ".fasta"
                file_paths["reference_genome_stat_path"] =  reference_dir_path + "/" + species.replace(" ","_") + "_genomeStats.txt"
                reference_genomes.append(file_paths["reference_genome_fasta_path"])

                #fetch the files
                response = urllib.request.urlopen(fasta_url)
                with gzip.GzipFile(fileobj=response) as remote_ref:
                    with open(file_paths["reference_genome_fasta_path"], 'wb') as local_ref:
                        shutil.copyfileobj(remote_ref, local_ref)
                urllib.request.urlretrieve(assembly_stat_url, file_paths["reference_genome_stat_path"])
                

    else: #throw an error if it contains contaminations
        print("Contaminated Genome assembly...resequencing required")
        raise Exception("contamination and mislabeling...crashing")
        
    #check to make sure we ONLY have ONE reference.
    if (len(reference_genomes) > 1 ):
        raise Exception ("there are multiple reference genomes")
    elif (len(reference_genomes) == 0):
        raise Exception ("no reference genome identified")
    
    #now we estimate our coverage using total reads and expected genome size
    expected_genome_size = result_parsers.parse_reference_genome_stats(file_paths["reference_genome_stat_path"])     #find expected genome size
    total_bp = result_parsers.parse_total_bp(file_paths["totalbp_path"])    #find total base count
    
    #calculate coverage
    coverage = total_bp / expected_genome_size
            
    if (coverage >= int(qc_cutoffs["coverage_cutoff"])):
        qc_verdicts["acceptable_coverage"] = True

    #time to assemble the reads then QC the assemblies
    print("step 2: genome assembly and QC")

    print("running shovill assembler")
    cmd = [script_path + "/job_scripts/shovill.sh",
           "--R1", R1, "--R2", R2,
           "--mincov", "3", "--minlen", "500", 
           "--output_dir", "/".join([outputDir, ID, "assembly"])]
    _ = execute(cmd, curDir)
    
    print("running busco")
    cmd = [script_path + "/job_scripts/busco.sh",
           "--input", "/".join([outputDir, ID, "assembly", "contigs.fa"]),
           "--database", buscodb, 
           "--output_dir", file_paths['busco_path']]
    _ = execute(cmd, curDir)

    print("running quast")
    cmd = [script_path + "/job_scripts/quast.sh",
           "--input", "/".join([outputDir, ID, "assembly", "contigs.fa"]),
           "--reference_genome", reference_genomes[0], 
           "--output_dir", file_paths['quast_path']]
    _ = execute(cmd, curDir)
    
    
    print("Parsing assembly results")
    #populate the busco and quast result object
    busco_results = result_parsers.parse_busco_result(file_paths["busco_path"] + "/short_summary_" + ID + ".busco.txt")
    quast_results = result_parsers.parse_quast_result(file_paths["quast_path"] + "/report.txt")    

    qc_verdicts["acceptable_busco_assembly_metrics"] = busco_results_qc_check(busco_results)
    qc_verdicts["acceptable_quast_assembly_metrics"] = quast_results_qc_check(quast_results)

    #print QC results to screen
    print("total bases: " + str(total_bp))
    print("expected genome size: " + str(expected_genome_size))
    print("coverage: " + str(coverage))
    print("")
    for key,value in qc_verdicts.items():
        print (str(key) + ": " + str(value))
        
    '''
    #print QC results into a txt
    print("Formatting the QC results")

    output.append("\n\n~~~~~~~QC summary~~~~~~~")
    output.append("Expected genome size: " + str(expected_genome_size))
    output.append("Estimated coverage: " + str(coverage))
    output.append("Expected isolate species: " + expectedSpecies)

    output.append("\nFastQC summary:")
    output.append("\nforward read qc:")
    for key, value in fastqcR1.items():
        output.append(key + ": " + value)
        if (value == "WARN" or value == "FAIL"):
            notes.append("FastQC: Forward read, " + key + " " + value)

    output.append("\nreverse read qc:")
    for key, value in fastqcR2.items():
        output.append(key + ": " + value)
        if (value == "WARN" or value == "FAIL"):
            notes.append("FastQC: Reverse read, " + key + " " + value)

    output.append("\nmash predicted genomes")
    for mash_hit in filtered_mash_hits:
        output.append(mash_hit['query_comment'])

    output.append("\nmash predicted plasmids")
    for mash_plasmid_hit in mash_plasmid_hits:
        output.append(mash_plasmid_hit['query_comment'])
    

    output.append("\nDetailed mash genome hits: ")
    for mash_hit in mash_hits:
        output.append(
            str(mash_hit['identity']) + '\t' +
            mash_hit['shared_hashes'] + '\t' +
            str(mash_hit['median_multiplicity']) + '\t' +
            str(mash_hit['p_value']) + '\t' +
            mash_hit['query_id'] + '\t' +
            mash_hit['query_comment']
        )
    
    output.append("\nDetailed mash plasmid hits: ")
    for mash_plasmid_hit in mash_plasmid_hits:
        output.append(
            str(mash_plasmid_hit['identity']) + '\t' +
            mash_plasmid_hit['shared_hashes'] + '\t' +
            str(mash_plasmid_hit['median_multiplicity']) + '\t' +
            str(mash_plasmid_hit['p_value']) + '\t' +
            mash_plasmid_hit['query_id'] + '\t' +
            mash_plasmid_hit['query_comment']
        )
    '''
    
if __name__ == "__main__":
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

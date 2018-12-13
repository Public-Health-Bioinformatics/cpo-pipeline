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

def httpGetFile(url, filepath=""):
    if (filepath == ""):
        return urllib.request.urlretrieve(url)
    else:
        urllib.request.urlretrieve(url, filepath)
        return True

def gunzip(file_path):
    with gzip.open(file_path, 'rb') as f:
        gzContent = f.read()
    with open(file_path.replace(".gz",""), 'wb') as out:
        out.write(gzContent)    
    os.remove(file_path)
    if (os.path.exists(file_path.replace(".gz",""))):
        return True
    else:
        raise Exception("Reference genome downloaded, but cannot be gunzipped")

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
    
    #parallelization, useless, these are hard coded to 8cores/64G RAM
    #parser.add_option("-t", "--threads", dest="threads", default=8, type="int", help="number of cpu to use")
    #parser.add_option("-p", "--memory", dest="memory", default=64, type="int", help="memory to use in GB")

    (options,args) = parser.parse_args()

    curDir = os.getcwd()
    outputDir = options.output
    mashdb = options.mashGenomeRefDB
    mashplasmiddb=options.mashPlasmidRefDB
    kraken2db = options.kraken2GenomeRefDB
    kraken2plasmiddb=options.kraken2PlasmidRefDB
    expectedSpecies = options.expectedSpecies
    #threads = options.threads
    #memory = options.memory
    tempDir = outputDir + "/shovillTemp"
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
        "mash_genome_path": (outputDir + "/qcResult/" + ID + "/" + "mashscreen.genome.tsv"),
        "mash_plasmid_path": (outputDir + "/qcResult/" + ID + "/" + "mashscreen.plasmid.tsv"),
        "fastqc_output_path":(outputDir + "/qcResult/" + ID + "/" + "fastqc"),
        "totalbp_path": (outputDir + "/qcResult/" + ID + "/" + "totalbp"),
        "reference_genome_fasta_path": (""), #built later
        "reference_genome_stat_path": (""), #built later
        "busco_path": (outputDir + "/assembly_qc/" + ID + "/" + ID + ".busco" + "/short_summary_" + ID + ".busco.txt"),
        "quast_path": (outputDir + "/assembly_qc/" + ID + "/" + ID + ".quast" + "/report.txt")
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
            qID = mash_hit['query_id'] #hit genome within mash results
            species_name_start = int(mash_hit['query_comment'].index(".")) + 3 #find the start of species name within query_comment column
            species_name_stop = int (mash_hit['query_comment'].index(",")) #find the end of the species name within query_comment column
            if (mash_hit['query_comment'].find("phiX") > -1):
                species = "PhiX" #phix
            else:
                species = str(mash_hit['query_comment'])[species_name_start: species_name_stop] #assign proper species name for reference genome file name
                # find gcf accession
                # TODO: document this or clean it up to be more readable
                gcf = (qID[:qID.find("_",5)]).replace("_","") #find the full gcf accession for ncbi FTP
                gcf = [gcf[i:i+3] for i in range(0, len(gcf), 3)] #break the gcf accession into k=3

                assembly = qID[:qID.find("_genomic.fna.gz")] #find the assembly name

                #build the urls
                fasta_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + gcf[0] + "/" + gcf[1] + "/" + gcf[2] + "/" + gcf[3] + "/" + assembly + "/" + qID #url to fasta
                assembly_stat_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + gcf[0] + "/" + gcf[1] + "/" + gcf[2] + "/" + gcf[3] + "/" + assembly + "/" + assembly + "_assembly_stats.txt" #url to assembly stat
                
                #build the save paths
                referencePath = os.path.abspath(outputDir + "/qcResult/" + ID + "/" + species.replace(" ",""))
                file_paths["reference_genome_fasta_path"] = (referencePath + ".fasta")
                file_paths["reference_genome_stat_path"] =  (referencePath + "_genomeStats.txt")
                reference_genomes.append(str(file_paths["reference_genome_fasta_path"]))

                #fetch the files
                httpGetFile(fasta_url, file_paths["reference_genome_fasta_path"] + ".gz") #fetch the fasta gz
                httpGetFile(assembly_stat_url, file_paths["reference_genome_stat_path"]) # fetch the genome stat
                
                #unzip the files
                gunzip(file_paths["reference_genome_fasta_path"] + ".gz")

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

    #run the assembly shell script.
    #input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome, 7=buscoDB
    cmd = [script_path + "/pipeline_assembly.sh", ID, R1, R2, outputDir, tempDir, reference_genomes[0], buscodb]
    result = execute(cmd, curDir)
    
    print("Parsing assembly results")
    #populate the busco and quast result object
    buscoResults = result_parsers.parse_busco_result(file_paths["busco_path"])
    quastResults = result_parsers.parse_quast_result(file_paths["quast_path"])
    
    #assembly QC logic    
    '''
    BUSCO PASS CRITERIA:
    1. complete singles > 90% of total genes
    2. complte duplicates < 90% of total genes
    '''
    if (float(buscoResults["complete_single"]) / float(buscoResults["total"]) >= float(qc_cutoffs["busco_complete_single_cutoff"]) and float(buscoResults["complete_duplicate"]) / float(buscoResults["total"]) <= float(qc_cutoffs["busco_complete_duplicate_cutoff"])):
        qc_verdicts["acceptable_busco_assembly_metrics"] = True

    '''
    QUAST PASS CRITERIA:
    1. total length vs reference length +-10%
    2. percent gc versus reference percent gc +- 5%
    3. genome fraction percent > 90
    '''   
    if ((float(quastResults["total_length"]) <= float(quastResults["reference_length"]) * (1 + float(qc_cutoffs["quast_assembly_length_cutoff"])) and float(quastResults["total_length"]) >= float(quastResults["reference_length"]) * (1 - float(qc_cutoffs["quast_assembly_length_cutoff"]))) and (float(quastResults["percent_GC"]) <= float(quastResults["reference_percent_GC"]) * (1+ float(qc_cutoffs["quast_percent_gc_cutoff"])) and float(quastResults["percent_GC"]) >= float(quastResults["reference_percent_GC"]) * (1 - float(qc_cutoffs["quast_percent_gc_cutoff"]))) and (float(quastResults["genome_fraction_percent"]) >= int(qc_cutoffs["genome_fraction_percent_cutoff"]))): 
        qc_verdicts["acceptable_quast_assembly_metrics"] = True

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

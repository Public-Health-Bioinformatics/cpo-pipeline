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
import pandas
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

def gunzip(inputpath="", outputpath=""):
    if (outputpath == ""):
        with gzip.open(inputpath, 'rb') as f:
            gzContent = f.read()
        return gzContent
    else:
        with gzip.open(inputpath, 'rb') as f:
            gzContent = f.read()
        with open(outputpath, 'wb') as out:
            out.write(gzContent)
        return True

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
    parser.add_option("-k", "--script-path", dest="scriptDir", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")

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
    scriptDir = options.scriptDir
    buscodb = options.buscodb
    ID = options.id
    R1 = options.R1
    R2 = options.R2

    
    notes = []
    #init the output list
    output = []
    
    print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nR1: " + R1 + "\nR2: " + R2)
    output.append(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nR1: " + R1 + "\nR2: " + R2)

    print("step 1: preassembly QC")

    print("running pipeline_qc.sh")
    #input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=mashgenomerefdb, $6=mashplasmidrefdb, $7=kraken2db, $8=kraken2plasmiddb
    cmd = [scriptDir + "/pipeline_qc.sh", ID, R1, R2, outputDir, mashdb, mashplasmiddb, kraken2db, kraken2plasmiddb]
    result = execute(cmd, curDir)

    print("Parsing the QC results")
    #parse read stats
    pathToMashLog = outputDir + "/qcResult/" + ID + "/" + "mash.log"
    pathToTotalBP = outputDir + "/qcResult/" + ID + "/" + "totalbp"
    stats = result_parsers.parse_read_stats(pathToMashLog, pathToTotalBP)

    #parse genome mash results
    pathToMashGenomeScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.genome.tsv"
    mashHits = result_parsers.parse_mash_genome_result(pathToMashGenomeScreenTSV, stats['size'], stats['depth'])

    # parse plasmid mash
    pathToMashPlasmidScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.plasmid.tsv"
    mashPlasmidHits = result_parsers.parse_mash_plasmid_result(pathToMashPlasmidScreenTSV, stats['size'], stats['depth'])

    # parse fastqc
    pathToFastQCR1 = outputDir + "/qcResult/" + ID + "/" + R1[R1.find(os.path.basename(R1)):R1.find(".")] + "_fastqc/summary.txt"
    pathToFastQCR2 = outputDir + "/qcResult/" + ID + "/" + R2[R2.find(os.path.basename(R2)):R2.find(".")] + "_fastqc/summary.txt"
    fastqcR1 = result_parsers.parse_fastqc_result(pathToFastQCR1)
    fastqcR2 = result_parsers.parse_fastqc_result(pathToFastQCR2)
    fastqc = {}
    fastqc["R1"]=fastqcR1
    fastqc["R2"]=fastqcR2
     
    # parse kraken2 result
    pathToKrakenResult = outputDir + "/qcResult/" + ID + "/kraken2.genome.report"
    krakenGenomes = result_parsers.parse_kraken_result(pathToKrakenResult)

    print("Formatting the QC results")
    multiple = False
    correctSpecies = False
    correctReference = ""

    output.append("\n\n~~~~~~~QC summary~~~~~~~")
    output.append("Estimated genome size: " + str(stats['size']))
    output.append("Estimated coverage: " + str(stats['depth']))
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

    output.append("\nKraken2 predicted species (>1%): ")
    for key in krakenGenomes:
        output.append(krakenGenomes[key]['name'])
    output.append("\nmash predicted genomes (within 300 of highest score): ")
    for key in mashHits:
        output.append(mashHits[key]['species'])
    output.append("\nmash predicted plasmids (within 100 of highest score): ")
    for key in mashPlasmidHits:
        output.append(mashPlasmidHits[key]['query_comment'])
    
    output.append("\nDetailed kraken genome hits: ")
    for key in krakenGenomes:
        output.append(krakenGenomes[key]['row'])
    output.append("\nDetailed mash genome hits: ")
    for key in mashHits:
        output.append(mashHits[key]['row'])
    output.append("\nDetailed mash plasmid hits: ")
    for key in mashPlasmidHits:
        output.append(mashPlasmidHits[key]['row'])

    #qcsummary
    output.append("\n\nQC Information:")

    present = False
    if (len(krakenGenomes) > 1):
        output.append("!!!Kraken2 predicted multiple species, possible contamination?")
        notes.append("Kraken2: multiple species, possible contamination.")
        #multiple = True
    elif (len(krakenGenomes) == 1):
        multiple = False
        
    for key in  krakenGenomes:
        if (krakenGenomes[key]['name'] == expectedSpecies):
            present = True

    if present:
        output.append("The expected species is predicted by kraken 2")
        #correctSpecies = True
    else:
        output.append("!!!The expected species is NOT predicted by kraken2, contamination? mislabeling?")
        notes.append("Kraken2: Not expected species. Possible contamination or mislabeling")

    if (stats['depth'] < 30):
        output.append("!!!Coverage is lower than 30. Estimated depth: " + str(stats['depth']))

    if (len(mashHits) > 1):
        output.append("!!!MASH predicted multiple species, possible contamination?")
        multiple=True
    elif (len(mashHits) < 1):
        output.append("!!!MASH had no hits, this is an unknown species")
        notes.append("Mash: Unknown Species")
                
    present=False
    for key in mashHits:
        spp = str(mashHits[key]['species'])
        if (spp.find(expectedSpecies) > -1 ):
            present=True
    if present:
        output.append("The expected species is predicted by mash")
        correctSpecies = True
        #notes.append("The expected species is predicted by mash")
    else:
        output.append("!!!The expected species is NOT predicted by mash, poor resolution? contamination? mislabeling?")
        notes.append("Mash: Not expected species. Possible resolution issues, contamination or mislabeling")

    if (len(mashPlasmidHits) == 0):
        output.append("!!!no plasmids predicted")
        notes.append("Mash: no plasmid predicted")

    #hack: throw exception if this analysis should not proceed due to contamination and mislabelling
    if (multiple and not correctSpecies):
        out = open(outputDir + "/summary/" + ID +".err", 'a')
        out.write('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')
        #raise Exception('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')
    print("Downloading reference genomes")

    referenceGenomes = []
    for key in mashHits:
        qID = mashHits[key]['query_id']

        # find gcf accession
        # TODO: document this or clean it up to be more readable
        gcf = (qID[:qID.find("_",5)]).replace("_","")
        gcf = [gcf[i:i+3] for i in range(0, len(gcf), 3)]
        
        assembly = qID[:qID.find("_genomic.fna.gz")]
        url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + gcf[0] + "/" + gcf[1] + "/" + gcf[2] + "/" + gcf[3] + "/" + assembly + "/" + qID
        referencePath = outputDir + "/qcResult/" + ID + "/" + key.replace(" ","")
        referenceGenomes.append(referencePath)

        httpGetFile(url, referencePath + ".gz")
        with gzip.open(referencePath + ".gz", 'rb') as f:
            gzContent = f.read()
        with open(referencePath, 'wb') as out:
            out.write(gzContent)
        os.remove(referencePath + ".gz")

    print("step 2: genome assembly and QC")
    correctAssembly = ""

    #input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome, 7=buscoDB
    #if (len(mashHits) == 1):
    if (len(referenceGenomes) > 1):
        for item in referenceGenomes:
            if (item.find(expectedSpecies.replace(" ", "")) > -1): #found the right genome
                correctAssembly = os.path.basename(item)
    else:
        correctAssembly = os.path.basename(referenceGenomes[0])
    if (correctAssembly == ""):
        raise Exception("no reference genome...crashing")

    if (not multiple and correctSpecies):
        print("Noncontaminated Genome assembly...")            
        cmd = [scriptDir + "/pipeline_assembly.sh", ID, R1, R2, outputDir, tempDir, referenceGenomes[0], buscodb, correctAssembly]
        result = execute(cmd, curDir)
    elif (multiple and correctSpecies):
        #input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome (csv, no spaces)	, 7=buscodb
        print("Contaminated Genome assembly...")
        cmd = [scriptDir + "/pipeline_assembly_contaminant.sh", ID, R1, R2, outputDir, tempDir, ",".join(referenceGenomes), buscodb, correctAssembly]
        result = execute(cmd, curDir)
    elif (multiple and not correctSpecies):
        print("Contaminated Genome assembly...No Correct Species Either")
        raise Exception("contamination and mislabeling...crashing")
        #cmd = [scriptDir + "/pipeline_assembly_contaminant.sh", ID, R1, R2, outputDir, tempDir, ",".join(referenceGenomes), buscodb, correctAssembly]
        #result = execute(cmd, curDir)
    elif (not multiple and not correctSpecies):
        print("Noncontaminated Genome assembly...No Correct species though")
        cmd = [scriptDir + "/pipeline_assembly.sh", ID, R1, R2, outputDir, tempDir, referenceGenomes[0], buscodb, correctAssembly]
        result = execute(cmd, curDir)

    print("Parsing assembly results")
    #get the correct busco and quast result file to parse
    correctAssembly = ""
    buscoPath = "" 
    quastPath = ""
    if ((not multiple and correctSpecies) or (not multiple and not correctSpecies)): #only 1 reference genome
        buscoPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + ".busco" + "/short_summary_" + ID + ".busco.txt")
        quastPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + ".quast" + "/report.tsv")
        correctAssembly = ID
    elif(multiple and correctSpecies): #multiple reference genome, need to find the one we care about
        for item in referenceGenomes:
            if (item.find(expectedSpecies.replace(" ", "")) > -1): #found the right genome
                buscoPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + "." + os.path.basename(item) + ".busco" + "/short_summary_" + ID + "." + os.path.basename(item) + ".busco.txt")
                quastPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + "." + os.path.basename(item) + ".quast/runs_per_reference/" + os.path.basename(item) + "/report.tsv")
                correctAssembly = ID + "." + os.path.basename(item)
    elif(multiple and not correctSpecies):
        raise Exception('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')

    if (buscoPath == "" or quastPath == ""):
        raise Exception('theres no reference genome for this sample for whatever reason...')

    #populate the busco and quast result object
    buscoResults = result_parsers.parse_busco_result(buscoPath)
    quastResults = result_parsers.parse_busco_result(quastPath)

if __name__ == "__main__":
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

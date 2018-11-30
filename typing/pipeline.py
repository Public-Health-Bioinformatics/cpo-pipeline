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
import numpy
import configparser

from parsers import result_parsers
from parsers import input_parsers

def read(path):
    return [line.rstrip('\n') for line in open(path)]

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

def ToJson(dictObject, outputPath):
    #outDir = outputDir + '/summary/' + ID + ".json/"
    #if not (os.path.exists(outDir)):
        #os.makedirs(outDir)
    #with open(outputPath, 'w') as f:
      #json.dump([ob.__dict__ for ob in dictObject.values()], f, ensure_ascii=False)
    return ""

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
    # if len(args) != 8:
        # parser.error("incorrect number of arguments, all 7 is required")
    curDir = os.getcwd()
    ID = str(options.id).lstrip().rstrip()
    assembly = options.assembly
    expected_species = options.expected_species
    script_path = options.script_path
    cardPath = options.card_path
    abricate_datadir = options.abricate_datadir
    abricate_cpo_plasmid_db = options.abricate_cpo_plasmid_db
    mlst_scheme_map_file = options.mlst_scheme_map_file
    # plasmidfinder = str(options.plasmidfinder).lstrip().rstrip()
    outputDir = options.output

    notes = []
    #init the output list
    output = []
    jsonOutput = []

    print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nAssembly: " + assembly)
    output.append(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nAssembly: " + assembly)

    #region call the typing script
    print("running pipeline_typing.sh")
    #input parameters: 1=id, 2=assembly, 3=output, 4=cardPath, 5=abricate_datadir, 6=abricate_cpo_plasmid_db
    cmd = [script_path + "/pipeline_typing.sh", ID, assembly, outputDir, cardPath, abricate_datadir, abricate_cpo_plasmid_db]
    result = execute(cmd, curDir)
    #endregion

    #region parse the mlst results
    
    print("step 3: parsing mlst, plasmid, and amr results")
    
    print("identifying MLST")
    mlst_report = outputDir + "/typing/" + ID + "/" + ID + ".mlst/" + ID + ".mlst" 
    mlstHits = result_parsers.parse_mlst_result(mlst_report)
    ToJson(mlstHits, "mlst.json")
    # TODO: Check that there is only one MLST result in the report, and handle
    #       cases where the report is malformed.
    mlstHit = mlstHits[0]
    mlst_scheme_map = input_parsers.parse_scheme_species_map(mlst_scheme_map_file)
    mlst_species = "Undefined"
    for scheme in mlst_scheme_map:
        if 'species' in scheme and scheme['scheme_id'] == mlstHit['scheme_id']:
            mlst_species = scheme['species']
            
    #endregion

    #region parse mobsuite, resfinder and rgi results
    print("identifying plasmid contigs and amr genes")

    #parse mobsuite results
    mob_recon_contig_report_path = outputDir + "/typing/" + ID + "/" + ID + ".recon/" + "contig_report.txt" 
    mob_recon_contig_report = result_parsers.parse_mob_recon_contig_report(mob_recon_contig_report_path)

    mob_recon_aggregate_report_path = outputDir + "/typing/" + ID + "/" + ID + ".recon/" + "mobtyper_aggregate_report.txt" 
    mSuitePlasmids = result_parsers.parse_mobsuite_plasmids(mob_recon_aggregate_report_path)
    ToJson(mSuitePlasmids, "mobsuitePlasmids.json")
    mob_recon_aggregate_report = result_parsers.parse_mob_recon_mobtyper_aggregate_report(mob_recon_aggregate_report_path)
    

    def extract_contig_num(contig_id):
        """
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

    def record_plasmid_contigs(mob_recon_contig_report):
        plasmid_contigs = []
        likely_plasmid_contigs = []
        for contig_report_record in mob_recon_contig_report:
            contig_num = extract_contig_num(contig_report_record['contig_id'])
            if contig_num not in plasmid_contigs and contig_num not in likely_plasmid_contigs:
                if contig_report_record['rep_type']:
                    plasmid_contigs.append(contig_num)
                else:
                    likely_plasmid_contigs.append(contig_num)
        return (plasmid_contigs, likely_plasmid_contigs)

    def record_plasmid_origins(mob_recon_contig_report):
        origins = []
        for contig_report_record in mob_recon_contig_report:
            if contig_report_record['rep_type']:
                if contig_report_record['rep_type'] not in origins:
                    print("Inserting " + contig_report_record['rep_type'] + " into origins")
                    origins.append(contig_report_record['rep_type'])
        return origins

    plasmid_contigs, likely_plasmid_contigs = record_plasmid_contigs(mob_recon_contig_report)

    origins = record_plasmid_origins(mob_recon_contig_report)
    
    #parse resfinder AMR results
    abricate = outputDir + "/resistance/" + ID + "/" + ID + ".cp"
    rFinder = result_parsers.parse_resfinder_result(abricate, plasmid_contigs, likely_plasmid_contigs)
    ToJson(rFinder, "resfinder.json")

    rgi = outputDir + "/resistance/" + ID + "/" + ID + ".rgi.txt"
    rgiAMR = result_parsers.parse_rgi_result(rgi, plasmid_contigs, likely_plasmid_contigs)
    ToJson(rgiAMR, "rgi.json")

    carbapenamases = []
    amrGenes = []
    for keys in rFinder:
        carbapenamases.append(rFinder[keys]['short_gene'] + "(" + rFinder[keys]['source'] + ")")
    for keys in rgiAMR:
        if (rgiAMR[keys]['drug_class'].find("carbapenem") > -1):
            if (rgiAMR[keys]['best_hit_aro'] not in carbapenamases):
                carbapenamases.append(rgiAMR[keys]['best_hit_aro'] + "(" + rgiAMR[keys]['source'] + ")")
        else:
            if (rgiAMR[keys]['best_hit_aro'] not in amrGenes):
                amrGenes.append(rgiAMR[keys]['best_hit_aro'] + "(" + rgiAMR[keys]['source'] + ")")
    #endregion

    #region output parsed mlst information
    print("formatting mlst outputs")
    output.append("\n\n\n~~~~~~~MLST summary~~~~~~~")
    output.append("MLST determined species: " + mlst_species)
    output.append("\nMLST Details: ")
    
    output.append("\nMLST information: ")
    if (mlst_species == expected_species):
        output.append("MLST determined species is the same as expected species")
        #notes.append("MLST determined species is the same as expected species")
    else:
        output.append("!!!MLST determined species is NOT the same as expected species, contamination? mislabeling?")
        notes.append("MLST: Not expected species. Possible contamination or mislabeling")

    #endregion

    #region output the parsed plasmid/amr results
    output.append("\n\n\n~~~~~~~~Plasmids~~~~~~~~\n")
    
    output.append("predicted plasmid origins: ")
    output.append(";".join(origins))

    output.append("\ndefinitely plasmid contigs")
    output.append(";".join(plasmid_contigs))
    
    output.append("\nlikely plasmid contigs")
    output.append(";".join(likely_plasmid_contigs))

    output.append("\nmob-suite prediction details: ")
    for mob_recon_contig_report_record in mob_recon_contig_report:
        output.append('\t'.join([str(x) for x in mob_recon_contig_report_record.values()]))

    output.append("\n\n\n~~~~~~~~AMR Genes~~~~~~~~\n")
    output.append("predicted carbapenamase Genes: ")
    output.append(",".join(carbapenamases))
    output.append("other RGI AMR Genes: ")
    for key in rgiAMR:
        output.append(rgiAMR[key]['best_hit_aro'] + "(" + rgiAMR[key]['source'] + ")")

    output.append("\nDetails about the carbapenamase Genes: ")
    for key in rFinder:
        output.append(rFinder[key]['row'])
    output.append("\nDetails about the RGI AMR Genes: ")
    for key in rgiAMR:
        output.append(rgiAMR[key]['row'])

    #write summary to a file
    summaryDir = outputDir + "/summary/" + ID
    os.makedirs(summaryDir, exist_ok=True)
    out = open(summaryDir + "/summary.txt", 'w')
    for item in output:
        out.write("%s\n" % item)


    #TSV output
    tsvOut = []
    tsvOut.append("ID\tExpected Species\tMLST Species\tSequence Type\tMLST Scheme\tCarbapenem Resistance Genes\tOther AMR Genes\tTotal Plasmids\tPlasmids ID\tNum_Contigs\tPlasmid Length\tPlasmid RepType\tPlasmid Mobility\tNearest Reference\tDefinitely Plasmid Contigs\tLikely Plasmid Contigs")
    #start with ID
    temp = ""
    temp += (ID + "\t")
    temp += expected_species + "\t"

    #move into MLST
    temp += mlst_species + "\t"
    temp += str(mlstHit['sequence_type']) + "\t"
    temp += mlstHit['scheme_id'] + "\t"
    
    #now onto AMR genes
    temp += ";".join(carbapenamases) + "\t"
    temp += ";".join(amrGenes) + "\t"

    #lastly plasmids
    temp+= str(len(mSuitePlasmids)) + "\t"
    plasmidID = ""
    contigs = ""
    lengths = ""
    rep_type = ""
    mobility = ""
    neighbour = ""
    for keys in mSuitePlasmids:
        plasmidID += str(mSuitePlasmids[keys]['mash_neighbor_cluster']) + ";"
        contigs += str(mSuitePlasmids[keys]['num_contigs']) + ";"
        lengths += str(mSuitePlasmids[keys]['total_length']) + ";"
        rep_type += str(mSuitePlasmids[keys]['rep_types']) + ";"
        mobility += str(mSuitePlasmids[keys]['predicted_mobility']) + ";"
        neighbour += str(mSuitePlasmids[keys]['mash_nearest_neighbor']) + ";"
    temp += plasmidID + "\t" + contigs + "\t" + lengths + "\t" + rep_type + "\t" + mobility + "\t" + neighbour + "\t"
    temp += ";".join(plasmid_contigs) + "\t"
    temp += ";".join(likely_plasmid_contigs)
    tsvOut.append(temp)

    summaryDir = outputDir + "/summary/" + ID
    out = open(summaryDir + "/summary.tsv", 'w')
    for item in tsvOut:
        out.write("%s\n" % item)
    #endregion

if __name__ == "__main__":
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

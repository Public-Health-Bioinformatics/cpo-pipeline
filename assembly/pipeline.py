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

from parsers.parsers import *

#region result objects
#define some objects to store values from results
#//TODO this is not the proper way of get/set private object variables. every value has manually assigned defaults intead of specified in init(). Also, use property(def getVar, def setVar).
class starFinders(object):
    def __init__(self):
        self.file = ""
        self.sequence = ""
        self.start = 0
        self.end = 0
        self.gene = ""
        self.shortGene = ""
        self.coverage = ""
        self.coverage_map = ""
        self.gaps = ""
        self.pCoverage = 100.00
        self.pIdentity = 100.00
        self.database = ""
        self.accession = ""
        self.product = ""
        self.source = "chromosome"
        self.row = ""

class PlasFlowResult(object):
    def __init__(self):
        self.sequence = ""
        self.length = 0
        self.label = ""
        self.confidence = 0
        self.usefulRow = ""
        self.row = ""

class MlstResult(object):
    def __init__(self):
        self.file = ""
        self.speciesID = ""
        self.seqType = 0
        self.scheme = ""
        self.species = ""
        self.row=""

class MashResult(object):
    def __init__(self):
        self.size = 0.0
        self.depth = 0.0
        self.identity = 0.0
        self.sharedHashes = ""
        self.medianMultiplicity = 0
        self.pvalue = 0.0
        self.queryID= ""
        self.queryComment = ""
        self.species = ""
        self.row = ""
        self.accession = ""
        self.gcf=""
        self.assembly=""

    def toDict(self): #doesnt actually work
        return dict((name, getattr(self, name)) for name in dir(self) if not name.startswith('__')) 

class fastqcResult(object):
    def __init__(self):
        self.basic = ""
        self.perBaseSequenceQuality = ""
        self.perTileSequenceQuality = ""
        self.perSequenceQualityScores = ""
        self.perBaseSequenceContent = ""
        self.perSequenceGCContent = ""
        self.perBaseNContent = ""
        self.sequenceLengthDistribution = ""
        self.sequenceDuplicationLevels = ""
        self.overrepresentedSequences = ""
        self.adapterContent = ""
        self.fastqcHTMLblob = ""
    def toDict(self):
        dict = {}
        dict["Basic Statistics"] = self.basic
        dict["Per base sequence quality"] = self.perBaseSequenceQuality
        dict["Per tile sequence quality"] = self.perTileSequenceQuality
        dict["Per sequence quality scores"] = self.perSequenceQualityScores
        dict["Per base sequence content"] = self.perBaseSequenceContent
        dict["Per sequence GC content"] = self.perSequenceGCContent
        dict["Per base N content"] = self.perBaseNContent
        dict["Sequence Length Distribution"] = self.sequenceLengthDistribution
        dict["Sequence Duplication Levels"] = self.sequenceDuplicationLevels
        dict["Overrepresented sequences"] = self.overrepresentedSequences
        dict["Adapter Content"] = self.adapterContent
        return dict
    def dataToTSV():
        dict = toDict()
        s = ""
        for key in dict:
            s += dict[key] + "\t"
        return s
    def headerToTSV():
        dict = toDict()
        s = ""
        for key in dict:
            s += key + "\t"
        return s

class KrakenResult(object):
    def __init__(self):
        self.fragPercent = -1.00
        self.fragCountRoot = 0
        self.fragCountTaxon = 0
        self.rankCode = ""
        self.taxID = -1
        self.name = ""
        self.row = ""

class buscoResult(object):
    def __init__(self):
        self.completeSingle = 0
        self.completeDuplicate = 0
        self.fragmented = 0
        self.mising = 0
        self.total = 781
        self.buscoResult = []

class quastResult(object):
    def __init__(self):
        self.contigsCount = 0
        self.largestContig = 0
        self.N50 = 0
        self.NG50 = 0
        self.L50 = 0
        self.LG50 = 0
        self.N75 = 0
        self.NG75 = 0
        self.L75 = 0
        self.LG75 = 0
        self.totalLength = 0
        self.referenceLength = 0
        self.GC = 0.00
        self.referenceGC = 0.00
        self.genomeFraction = 0.00
        self.duplicationRatio = 0.00
        self.data = []

class mobsuiteResult(object):
    def __init__(self):
        self.file_id = ""
        self.cluster_id	= ""
        self.contig_id	= ""
        self.contig_num = 0
        self.contig_length	= 0
        self.circularity_status	= ""
        self.rep_type	= ""
        self.rep_type_accession = ""	
        self.relaxase_type	= ""
        self.relaxase_type_accession = ""	
        self.mash_nearest_neighbor	 = ""
        self.mash_neighbor_distance	= 0.00
        self.repetitive_dna_id	= ""
        self.match_type	= ""
        self.score	= 0
        self.contig_match_start	= 0
        self.contig_match_end = 0
        self.row = ""

class RGIResult(object):
    def __init__(self):
        self.ORF_ID	= ""
        self.Contig	= ""
        self.Start	= -1
        self.Stop	= -1
        self.Orientation = ""	
        self.Cut_Off	= ""
        self.Pass_Bitscore	= 100000
        self.Best_Hit_Bitscore	= 0.00
        self.Best_Hit_ARO	= ""
        self.Best_Identities	= 0.00
        self.ARO = 0
        self.Model_type	= ""
        self.SNPs_in_Best_Hit_ARO	= ""
        self.Other_SNPs	= ""
        self.Drug_Class	= ""
        self.Resistance_Mechanism	= ""
        self.AMR_Gene_Family	= ""
        self.Predicted_DNA	= ""
        self.Predicted_Protein	= ""
        self.CARD_Protein_Sequence	= ""
        self.Percentage_Length_of_Reference_Sequence	= 0.00
        self.ID	= ""
        self.Model_ID = 0
        self.source = ""
        self.row = ""

#endregion

#region useful functions
def read(path):
    return [line.rstrip('\n') for line in open(path)]

def execute(command):
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
    outDir = outputDir + '/summary/' + ID + ".json/"
    #if not (os.path.exists(outDir)):
        #os.makedirs(outDir)
    #with open(outDir + outputPath, 'w') as f:
        #json.dump([ob.__dict__ for ob in dictObject.values()], f, ensure_ascii=False)
#endregion

#region functions to parse result files
def ParseKrakenResult(pathToKrakenResult):
    _krakenGenomes = {}
    kraken = pandas.read_csv(pathToKrakenResult, delimiter='\t', header=None) #read  the kraken2report.tsv
    kraken = kraken.loc[(kraken[3] == "S")] #find all the rows with species level information and discard the rest
    for i in range(len(kraken.index)):
        if (float((str(kraken.iloc[i, 0])).strip("%")) > 1.00): #find all species above 1%
            kg = KrakenResult() 
            kg.fragPercent = float((str(kraken.iloc[i, 0])).strip("%"))
            kg.fragCountRoot = int(kraken.iloc[i,1])
            kg.fragCountTaxon = int (kraken.iloc[i,2])
            kg.rankCode = kraken.iloc[i,3]
            kg.taxID = kraken.iloc[i,4]
            kg.name = kraken.iloc[i,5].strip()
            kg.row = "\t".join(str(x) for x in kraken.iloc[i].tolist())
            _krakenGenomes[kg.name]=(kg) #throw the kraken result object in to a list
    return _krakenGenomes

def ParseFastQCResult(pathToR1qc, pathToR2qc, ID, R1, R2):
    fastqc = pandas.read_csv(pathToR1qc + "summary.txt", delimiter='\t', header=None)
    fastqc = fastqc[0].tolist()

    _fastqcR1 = fastqcResult()
    _fastqcR1.basic = fastqc[0]
    _fastqcR1.perBaseSequenceQuality = fastqc[1]
    _fastqcR1.perTileSequenceQuality = fastqc[2]
    _fastqcR1.perSequenceQualityScores = fastqc[3]
    _fastqcR1.perBaseSequenceContent = fastqc[4]
    _fastqcR1.perSequenceGCContent = fastqc[5]
    _fastqcR1.perBaseNContent = fastqc[6]
    _fastqcR1.sequenceLengthDistribution = fastqc[7]
    _fastqcR1.sequenceDuplicationLevels = fastqc[8]
    _fastqcR1.overrepresentedSequences = fastqc[9]
    _fastqcR1.adapterContent = fastqc[10]
    with open(pathToR1qc + "fastqc_report.html", "r") as input:
        _fastqcR1.fastqcHTMLblob = input.read()

    fastqc = pandas.read_csv(pathToR2qc+ "summary.txt", delimiter='\t', header=None)
    fastqc = fastqc[0].tolist()

    _fastqcR2 = fastqcResult()
    _fastqcR2.basic = fastqc[0]
    _fastqcR2.perBaseSequenceQuality = fastqc[1]
    _fastqcR2.perTileSequenceQuality = fastqc[2]
    _fastqcR2.perSequenceQualityScores = fastqc[3]
    _fastqcR2.perBaseSequenceContent = fastqc[4]
    _fastqcR2.perSequenceGCContent = fastqc[5]
    _fastqcR2.perBaseNContent = fastqc[6]
    _fastqcR2.sequenceLengthDistribution = fastqc[7]
    _fastqcR2.sequenceDuplicationLevels = fastqc[8]
    _fastqcR2.overrepresentedSequences = fastqc[9]
    _fastqcR2.adapterContent = fastqc[10]
    with open(pathToR2qc + "fastqc_report.html", "r") as input:
        _fastqcR2.fastqcHTMLblob = input.read()

    return _fastqcR1, _fastqcR2

def ParseMashGenomeResult(pathToMashScreen, size, depth):
    _mashHits = {}
    _PhiX = False
    
    mashScreen = pandas.read_csv(pathToMashScreen, delimiter='\t', header=None) #read mash report
    #0.998697        973/1000        71      0       GCF_000958965.1_matepair4_genomic.fna.gz        [59 seqs] NZ_LAFU01000001.1 Klebsiella pneumoniae strain CDPH5262 contig000001, whole genome shotgun sequence [...]
    #parse mash result, using winner takes all
    scores = mashScreen[1].values 
    scoreCutoff = int(scores[0][:scores[0].index("/")]) - 300 #find the cut off value to use for filtering the report (300 below max)
    index = 0
    #find hits with score within top 300
    for score in scores:
        parsedScore = int(score[:score.index("/")])
        if parsedScore >= scoreCutoff:
            index+=1
        else:
            break

    #parse what the species are.
    for i in range(index):
        mr = MashResult()
        mr.size = size
        mr.depth = depth
        mr.identity = float(mashScreen.ix[i, 0])
        mr.sharedHashes = mashScreen.ix[i, 1]
        mr.medianMultiplicity = int(mashScreen.ix[i, 2])
        mr.pvalue = float(mashScreen.ix[i, 3])
        mr.queryID = mashScreen.ix[i, 4]
        mr.queryComment = mashScreen.ix[i, 5]
        mr.accession = mr.queryComment
        mr.row = "\t".join(str(x) for x in mashScreen.ix[i].tolist())
        qID = mr.queryID
        #find gcf accession
        gcf = (qID[:qID.find("_",5)]).replace("_","")
        gcf = [gcf[i:i+3] for i in range(0, len(gcf), 3)]
        mr.gcf = gcf 
        #find assembly name
        mr.assembly = qID[:qID.find("_genomic.fna.gz")]
        #score = mashScreen.iloc[[i]][1]

        if (mr.queryComment.find("phiX") > -1): #theres phix in top hits, just ignore
            _PhiX=True
            mr.species = "PhiX"
        else: #add the non-phix hits to a list
            start = int(mr.queryComment.index(".")) + 3
            stop = int (mr.queryComment.index(","))
            mr.species = str(mr.queryComment)[start: stop]
            _mashHits[mr.species]=mr
    return _mashHits, _PhiX

def ParseMashPlasmidResult(pathToMashScreen, size, depth):
    mashScreen = pandas.read_csv(pathToMashScreen, delimiter='\t', header=None)

    #parse mash result, using winner takes all
    scores = mashScreen[1].values 
    scoreCutoff = int(scores[0][:scores[0].index("/")]) - 100
    index = 0

    #find hits with score >850
    for score in scores:
        parsedScore = int(score[:score.index("/")])
        if parsedScore >= scoreCutoff:
            index+=1
        else:
            break

    _mashPlasmidHits = {}
    #parse what the species are.
    for i in range(index):
        mr = MashResult()
        mr.size = size
        mr.depth = depth
        mr.identity = float(mashScreen.ix[i, 0])
        mr.sharedHashes = mashScreen.ix[i, 1]
        mr.medianMultiplicity = int(mashScreen.ix[i, 2])
        mr.pvalue = float(mashScreen.ix[i, 3])
        mr.queryID = mashScreen.ix[i, 4] #accession
        mr.queryComment = mashScreen.ix[i, 5] #what is it
        mr.accession = mr.queryID[4:len(mr.queryID)-1]
        mr.row = "\t".join(str(x) for x in mashScreen.ix[i].tolist())
        #score = mashScreen.iloc[[i]][1]
        mr.species = ""
        if (mr.identity >= 0.97):
            _mashPlasmidHits[mr.accession] = mr
    return _mashPlasmidHits

def ParseReadStats(pathToMashLog, pathToTotalBp):
    for s in read(pathToMashLog):
        if (s.find("Estimated genome size:") > -1 ):
            _size = float(s[s.index(": ")+2:])
    _totalbp = float(read(pathToTotalBp)[0])
    _depth = _totalbp / _size
    _depth = float(format(_depth, '.2f'))
    return _size,_depth

def ParseBuscoResult(pathToBuscoResult):
    buscoOutput = read(pathToBuscoResult)
    if (len(buscoOutput) > 0):
        br = buscoResult()
        br.completeSingle = int(buscoOutput[10].split("\t")[1])
        br.completeDuplicate = int(buscoOutput[11].split("\t")[1])
        br.fragmented = int(buscoOutput[12].split("\t")[1])
        br.mising = int(buscoOutput[13].split("\t")[1])
        br.total = int(buscoOutput[14].split("\t")[1])
        br.buscoResult = buscoOutput
        _buscoResults = br
    else: #it should never be <1.... i think
        raise Exception("x011 this shouldnt happen...i think")
    return _buscoResults

def ParseQuastResult(pathToQuastResult):
    qResult = read(pathToQuastResult)
    qr = quastResult()

    qr.contigsCount = int(qResult[int([i for i,x in enumerate(qResult) if x.find("# contigs\t") > -1][-1])].split("\t")[1])
    qr.largestContig = int(qResult[int([i for i,x in enumerate(qResult) if x.find("Largest contig\t") > -1][-1])].split("\t")[1])
    qr.N50 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("N50\t") > -1][-1])].split("\t")[1])
    qr.NG50 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("NG50\t") > -1][-1])].split("\t")[1])
    qr.L50 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("L50\t") > -1][-1])].split("\t")[1])
    qr.LG50 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("LG50\t") > -1][-1])].split("\t")[1])
    qr.N75 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("N75\t") > -1][-1])].split("\t")[1])
    qr.NG75 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("NG75\t") > -1][-1])].split("\t")[1])
    qr.L75 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("L75\t") > -1][-1])].split("\t")[1])
    qr.LG75 = int(qResult[int([i for i,x in enumerate(qResult) if x.find("LG75\t") > -1][-1])].split("\t")[1])
    qr.totalLength = int(qResult[int([i for i,x in enumerate(qResult) if x.find("Total length\t") > -1][-1])].split("\t")[1])
    qr.referenceLength = int(qResult[int([i for i,x in enumerate(qResult) if x.find("Reference length\t") > -1][-1])].split("\t")[1])
    qr.GC = float(qResult[int([i for i,x in enumerate(qResult) if x.find("GC (%)\t") > -1][-1])].split("\t")[1])
    qr.referenceGC = float(qResult[int([i for i,x in enumerate(qResult) if x.find("Reference GC (%)\t") > -1][-1])].split("\t")[1])
    qr.genomeFraction = float(qResult[int([i for i,x in enumerate(qResult) if x.find("Genome fraction (%)\t") > -1][-1])].split("\t")[1])
    qr.duplicationRatio = float(qResult[int([i for i,x in enumerate(qResult) if x.find("Duplication ratio\t") > -1][-1])].split("\t")[1])
    qr.data = qResult

    return qr
#endregion

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
    #if len(args) != 8:
        #parser.error("incorrect number of arguments, all 7 is required")

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

    #region step 1 pre assembly qc
    print("step 1: preassembly QC")

    #region call the qc script
    print("running pipeline_qc.sh")
    #input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=mashgenomerefdb, $6=mashplasmidrefdb, $7=kraken2db, $8=kraken2plasmiddb
    cmd = [scriptDir + "/pipeline_qc.sh", ID, R1, R2, outputDir, mashdb, mashplasmiddb, kraken2db, kraken2plasmiddb]
    result = execute(cmd)
    #endregion

    print("Parsing the QC results")
    #region parse results
    #parse read stats
    pathToMashLog = outputDir + "/qcResult/" + ID + "/" + "mash.log"
    pathToTotalBP = outputDir + "/qcResult/" + ID + "/" + "totalbp"
    size, depth = ParseReadStats(pathToMashLog, pathToTotalBP)
    stats = {}
    stats["size"] = size
    stats["depth"] = depth
    ToJson(stats, "readStats.json")

    #parse genome mash results
    pathToMashGenomeScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.genome.tsv"
    mashHits, PhiX = ParseMashGenomeResult(pathToMashGenomeScreenTSV, size, depth)
    ToJson(mashHits, "mashGenomeHit.json")

    # parse plasmid mash
    pathToMashPlasmidScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.plasmid.tsv"
    mashPlasmidHits = ParseMashPlasmidResult(pathToMashPlasmidScreenTSV, size, depth)
    ToJson(mashPlasmidHits, "mashPlasmidHits.json")

    # parse fastqc
    pathToFastQCR1 = outputDir + "/qcResult/" + ID + "/" + R1[R1.find(os.path.basename(R1)):R1.find(".")] + "_fastqc/"
    pathToFastQCR2 = outputDir + "/qcResult/" + ID + "/" + R2[R2.find(os.path.basename(R2)):R2.find(".")] + "_fastqc/"
    fastqcR1,fastqcR2 = ParseFastQCResult(pathToFastQCR1, pathToFastQCR2, ID, R1, R2)
    fastqc = {}
    fastqc["R1"]=fastqcR1
    fastqc["R2"]=fastqcR1
    ToJson(fastqc, "fastqc.json")
     
    # parse kraken2 result
    pathToKrakenResult = outputDir + "/qcResult/" + ID + "/kraken2.genome.report"
    krakenGenomes = ParseKrakenResult(pathToKrakenResult)
    ToJson(krakenGenomes, "krakenGenomeHits.json")
    #pathToKrakenPlasmidResult = outputDir + "/qcResult/" + ID + "/kraken2.plasmid.report"
    
    #endregion

    print("Formatting the QC results")
    #region output parsed results
    multiple = False
    correctSpecies = False
    correctReference = ""

    output.append("\n\n~~~~~~~QC summary~~~~~~~")
    output.append("Estimated genome size: " + str(size))
    output.append("Estimated coverage: " + str(depth))
    output.append("Expected isolate species: " + expectedSpecies)

    output.append("\nFastQC summary:")
    output.append("\nforward read qc:")
    for key in fastqcR1.toDict():
        output.append(key + ": " + fastqcR1.toDict()[key])
        if (fastqcR1.toDict()[key] == "WARN" or fastqcR1.toDict()[key] == "FAIL"):
            notes.append("FastQC: Forward read, " + key + " " + fastqcR1.toDict()[key])
    output.append("\nreverse read qc:")
    for key in fastqcR2.toDict():
        output.append(key + ": " + fastqcR2.toDict()[key])
        if (fastqcR2.toDict()[key] == "WARN" or fastqcR2.toDict()[key] == "FAIL"):
            notes.append("FastQC: Reverse read, " + key + " " + fastqcR2.toDict()[key])

    output.append("\nKraken2 predicted species (>1%): ")
    for key in krakenGenomes:
        output.append(krakenGenomes[key].name)
    output.append("\nmash predicted genomes (within 300 of highest score): ")
    for key in mashHits:
        output.append(mashHits[key].species)
    output.append("\nmash predicted plasmids (within 100 of highest score): ")
    for key in mashPlasmidHits:
        output.append(mashPlasmidHits[key].queryComment)
    
    output.append("\nDetailed kraken genome hits: ")
    for key in krakenGenomes:
        output.append(krakenGenomes[key].row)
    output.append("\nDetailed mash genome hits: ")
    for key in mashHits:
        output.append(mashHits[key].row)
    output.append("\nDetailed mash plasmid hits: ")
    for key in mashPlasmidHits:
        output.append(mashPlasmidHits[key].row)

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
        if (krakenGenomes[key].name == expectedSpecies):
            present = True

    if present:
        output.append("The expected species is predicted by kraken 2")
        #correctSpecies = True
    else:
        output.append("!!!The expected species is NOT predicted by kraken2, contamination? mislabeling?")
        notes.append("Kraken2: Not expected species. Possible contamination or mislabeling")

    if (depth < 30):
        output.append("!!!Coverage is lower than 30. Estimated depth: " + str(depth))

    if(PhiX):
        output.append("!!!PhiX contamination, probably nothing to worry about")
        notes.append("PhiX contamination")

    if (len(mashHits) > 1):
        output.append("!!!MASH predicted multiple species, possible contamination?")
        multiple=True
    elif (len(mashHits) < 1):
        output.append("!!!MASH had no hits, this is an unknown species")
        notes.append("Mash: Unknown Species")
                
    present=False
    for key in mashHits:
        spp = str(mashHits[key].species)
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
    #endregion

    #hack: throw exception if this analysis should not proceed due to contamination and mislabelling
    if (multiple and not correctSpecies):
        out = open(outputDir + "/summary/" + ID +".err", 'a')
        out.write('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')
        #raise Exception('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')
    print("Downloading reference genomes")
    #region identify and download a reference genome with mash
    referenceGenomes = []
    for key in mashHits:
        qID = mashHits[key].queryID
        gcf = mashHits[key].gcf
        assembly = mashHits[key].assembly
        url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + gcf[0] + "/" + gcf[1] + "/" + gcf[2] + "/" + gcf[3] + "/" + assembly + "/" + qID
        referencePath = outputDir + "/qcResult/" + ID + "/" + mashHits[key].species.replace(" ","")
        referenceGenomes.append(referencePath)

        httpGetFile(url, referencePath + ".gz")
        with gzip.open(referencePath + ".gz", 'rb') as f:
            gzContent = f.read()
        with open(referencePath, 'wb') as out:
            out.write(gzContent)
        os.remove(referencePath + ".gz")
    
    #endregion

    #endregion

    #region step2, assembly and qc
    print("step 2: genome assembly and QC")

    #region call the script
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
        result = execute(cmd)
    elif (multiple and correctSpecies):
        #input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome (csv, no spaces)	, 7=buscodb
        print("Contaminated Genome assembly...")
        cmd = [scriptDir + "/pipeline_assembly_contaminant.sh", ID, R1, R2, outputDir, tempDir, ",".join(referenceGenomes), buscodb, correctAssembly]
        result = execute(cmd)
    elif (multiple and not correctSpecies):
        print("Contaminated Genome assembly...No Correct Species Either")
        raise Exception("contamination and mislabeling...crashing")
        #cmd = [scriptDir + "/pipeline_assembly_contaminant.sh", ID, R1, R2, outputDir, tempDir, ",".join(referenceGenomes), buscodb, correctAssembly]
        #result = execute(cmd)
    elif (not multiple and not correctSpecies):
        print("Noncontaminated Genome assembly...No Correct species though")
        cmd = [scriptDir + "/pipeline_assembly.sh", ID, R1, R2, outputDir, tempDir, referenceGenomes[0], buscodb, correctAssembly]
        result = execute(cmd)
    #endregion

    print("Parsing assembly results")
    #region parse assembly_qc from quast and busco
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
    elif(multiple and not correctSpecies): #the most fked up case, dont even bother? just crash the damn program for now
        raise Exception('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')

    if (buscoPath == "" or quastPath == ""):
        raise Exception('theres no reference genome for this sample for whatever reason...')

    #populate the busco and quast result object
    buscoResults = ParseBuscoResult(buscoPath)
    quastResults = ParseQuastResult(quastPath)

    #endregion

    #endregion


if __name__ == "__main__":
    start = time.time()
    print("Starting workflow...")
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end - start) + " seconds")

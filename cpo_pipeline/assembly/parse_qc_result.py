#!/usr/bin/env python

import pandas
import json

def ParseReadStats(pathToMashLog, pathToTotalBp):
    for s in read(pathToMashLog):
        if (s.find("Estimated genome size:") > -1 ):
            _size = float(s[s.index(": ")+2:])
    _totalbp = float(read(pathToTotalBp)[0])
    _depth = _totalbp / _size
    _depth = float(format(_depth, '.2f'))
    return _size,_depth

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
    
def ParseMashGenomeResult(pathToMashScreen, size, depth):
    _mashHits = {} #***********************
    _PhiX = False #***********************
    
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

    _mashPlasmidHits = {} #***********************
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

class KrakenResult(object):
    def __init__(self):
        self.fragPercent = -1.00
        self.fragCountRoot = 0
        self.fragCountTaxon = 0
        self.rankCode = ""
        self.taxID = -1
        self.name = ""
        self.row = ""

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

def ToJson(dictObject, outputPath):
    outDir = outputDir + '/summary/' + ID + ".json/"
    #if not (os.path.exists(outDir)):
        #os.makedirs(outDir)
    #with open(outDir + outputPath, 'w') as f:
        #json.dump([ob.__dict__ for ob in dictObject.values()], f, ensure_ascii=False)

print("Parsing the QC results")
#region parse results
#parse read stats
pathToMashLog = outputDir + "/qcResult/" + ID + "/" + "mash.log"
pathToTotalBP = outputDir + "/qcResult/" + ID + "/" + "totalbp"
size, depth = ParseReadStats(pathToMashLog, pathToTotalBP) #**************
stats = {}
stats["size"] = size
stats["depth"] = depth
ToJson(stats, "readStats.json")

#parse genome mash results
pathToMashGenomeScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.genome.tsv"
mashHits, PhiX = ParseMashGenomeResult(pathToMashGenomeScreenTSV, size, depth) #***********************
ToJson(mashHits, "mashGenomeHit.json")

# parse plasmid mash
pathToMashPlasmidScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.plasmid.tsv"
mashPlasmidHits = ParseMashPlasmidResult(pathToMashPlasmidScreenTSV, size, depth) #**************
ToJson(mashPlasmidHits, "mashPlasmidHits.json")

# parse fastqc
pathToFastQCR1 = outputDir + "/qcResult/" + ID + "/" + R1[R1.find(os.path.basename(R1)):R1.find(".")] + "_fastqc/"
pathToFastQCR2 = outputDir + "/qcResult/" + ID + "/" + R2[R2.find(os.path.basename(R2)):R2.find(".")] + "_fastqc/"
fastqcR1,fastqcR2 = ParseFastQCResult(pathToFastQCR1, pathToFastQCR2, ID, R1, R2) #**************
fastqc = {}
fastqc["R1"]=fastqcR1
fastqc["R2"]=fastqcR1
ToJson(fastqc, "fastqc.json")
     
# parse kraken2 result
pathToKrakenResult = outputDir + "/qcResult/" + ID + "/kraken2.genome.report"
krakenGenomes = ParseKrakenResult(pathToKrakenResult) #**************
ToJson(krakenGenomes, "krakenGenomeHits.json")
#pathToKrakenPlasmidResult = outputDir + "/qcResult/" + ID + "/kraken2.plasmid.report"

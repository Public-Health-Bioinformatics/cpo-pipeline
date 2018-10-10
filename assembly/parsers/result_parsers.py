def parse_kraken_result(pathToKrakenResult):
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

def parse_fastqc_result(pathToR1qc, pathToR2qc, ID, R1, R2):
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

def parse_mash_genome_result(pathToMashScreen, size, depth):
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

def parse_mash_plasmid_result(pathToMashScreen, size, depth):
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

def parse_read_stats(pathToMashLog, pathToTotalBp):
    for s in read(pathToMashLog):
        if (s.find("Estimated genome size:") > -1 ):
            _size = float(s[s.index(": ")+2:])
    _totalbp = float(read(pathToTotalBp)[0])
    _depth = _totalbp / _size
    _depth = float(format(_depth, '.2f'))
    return _size,_depth

def parse_busco_result(pathToBuscoResult):
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

def parse_quast_result(pathToQuastResult):
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


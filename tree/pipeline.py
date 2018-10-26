#!/usr/bin/env python

#$ -S /usr/bin/env python
#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 1      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log
#$ -m ea


# This script is a wrapper for module two, part 2: tree rendering of the cpo_workflow to render dendrograms.
# It uses snippy for core genome SNV calling and alignment, clustalw to generate a NJ tree and ete3 to render the dendrogram

#   >python cpo_galaxy_tree.py -t /path/to/tree.ph -d /path/to/distance/matrix -m /path/to/metadata

#	<requirements>
#		<requirement type="package" version="0.23.4">pandas</requirement>
#		<requirement type="package" version="3.6">python</requirement>
#       <requirement type="package" version="3.1.1">ete3</requirement>
#		<requirement type="package" version="5.6.0">pyqt</requirement>
#		<requirement type="package" version="5.6.2">qt</requirement>
#  </requirements>

import subprocess
import pandas #conda pandas
import optparse
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import datetime
import sys
import time
import urllib.request
import gzip
import collections
import json
import numpy #conda numpy
import ete3 as e #conda ete3 3.1.1**** >requires pyqt5

from parsers import result_parsers

#parses some parameters
parser = optparse.OptionParser("Usage: %prog [options] arg1 arg2 ...")
parser.add_option("-t", "--tree", dest="treePath", type="string", default="./pipelineTest/tree.txt", help="absolute file path to phylip tree")    
parser.add_option("-d", "--distance", dest="distancePath", type="string", default="./pipelineTest/dist.tabular", help="absolute file path to distance matrix")
parser.add_option("-m", "--metadata", dest="metadataPath", type="string", default="./pipelineTest/metadata.tabular",help="absolute file path to metadata file")
parser.add_option("-o", "--output_file", dest="outputFile", type="string", default="tree.png", help="Output graphics file. Use ending 'png', 'pdf' or 'svg' to specify file format.")

# sensitive data adder
parser.add_option("-p", "--sensitive_data", dest="sensitivePath", type="string", default="", help="Spreadsheet (CSV) with sensitive metadata")
parser.add_option("-c", "--sensitive_cols", dest="sensitiveCols", type="string", default="", help="CSV list of column names from sensitive metadata spreadsheet to use as labels on dendrogram")
parser.add_option("-b", "--bcid_column", dest="bcidCol", type="string", default="BCID", help="Column name of BCID in sensitive metadata file")
parser.add_option("-n", "--missing_value", dest="naValue", type="string", default="NA", help="Value to write for missing data.")

(options,args) = parser.parse_args()
treePath = str(options.treePath).lstrip().rstrip()
distancePath = str(options.distancePath).lstrip().rstrip()
metadataPath = str(options.metadataPath).lstrip().rstrip()

sensitivePath = str(options.sensitivePath).lstrip().rstrip()
sensitiveCols = str(options.sensitiveCols).lstrip().rstrip()
outputFile = str(options.outputFile).lstrip().rstrip()
bcidCol = str( str(options.bcidCol).lstrip().rstrip() ) 
naValue = str( str(options.naValue).lstrip().rstrip() ) 


#region result objects
#define some objects to store values from results
#//TODO this is not the proper way of get/set private object variables. every value has manually assigned defaults intead of specified in init(). Also, use property(def getVar, def setVar).

class SensitiveMetadata(object):
    def __init__(self):
        x = pandas.read_csv( sensitivePath )
        col_names = [ s for s in sensitiveCols.split(',')] # convert to 0 offset
        if not bcidCol in col_names:
            col_names.append( bcidCol )
        all_cols = [ str(col) for col in x.columns ]
        col_idxs = [ all_cols.index(col) for col in col_names ]
        self.sensitive_data = x.iloc[:, col_idxs]
    def get_columns(self):
        cols = [ str(x) for x in self.sensitive_data.columns ]
        return cols
    def get_value( self, bcid, column_name ): # might be nice to get them all in single call via an input list of bcids ... for later
        bcids= list( self.sensitive_data.loc[:, bcidCol ] ) # get the list of all BCIDs  in sensitive metadata
        if not bcid in bcids:
            return naValue
        else:
            row_idx = bcids.index( bcid )                        # lookup the row for this BCID
        return self.sensitive_data.loc[ row_idx, column_name ]  # return the one value based on the column (col_idx) and this row



def read(path): #read in a text file to a list
    return [line.rstrip('\n') for line in open(path)]

def execute(command): #subprocess.popen call bash command
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

def httpGetFile(url, filepath=""): #download a file from the web
    if (filepath == ""):
        return urllib.request.urlretrieve(url)
    else:
        urllib.request.urlretrieve(url, filepath)
        return True
    
def gunzip(inputpath="", outputpath=""): #gunzip
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
    
def addFace(name): #function to add a facet to a tree
    #if its the reference branch, populate the faces with column headers
    face = e.faces.TextFace(name,fsize=10,tight_text=True)
    face.border.margin = 5
    face.margin_right = 5
    face.margin_left = 5
    return face


def main():
    if len(sensitivePath)>0:
        sensitive_meta_data = SensitiveMetadata()

    metadata = result_parsers.parse_workflow_results(metadataPath)
    distance = read(distancePath)
    treeFile = "".join(read(treePath))

    distanceDict = {} #store the distance matrix as rowname:list<string>
    for i in range(len(distance)):
        temp = distance[i].split("\t")
        distanceDict[temp[0]] = temp[1:]

    #region create box tree
    #region step5: tree construction
    treeFile = "".join(read(treePath))
    t = e.Tree(treeFile)
    t.set_outgroup(t&"Reference")

    #set the tree style
    ts = e.TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.scale = 2000 #pixel per branch length unit
    ts.branch_vertical_margin = 15 #pixel between branches
    style2 = e.NodeStyle()
    style2["fgcolor"] = "#000000"
    style2["shape"] = "circle"
    style2["vt_line_color"] = "#0000aa"
    style2["hz_line_color"] = "#0000aa"
    style2["vt_line_width"] = 2
    style2["hz_line_width"] = 2
    style2["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
    style2["hz_line_type"] = 0
    for n in t.traverse():
        n.set_style(style2)

    #find the plasmid origins
    plasmidIncs = {}
    for key in metadata:
        for plasmid in metadata[key].plasmids:
            for inc in plasmid.PlasmidRepType.split(","):
                if (inc.lower().find("inc") > -1):
                    if not (inc in plasmidIncs):
                        plasmidIncs[inc] = [metadata[key].ID]
                    else:
                        if metadata[key].ID not in plasmidIncs[inc]:
                            plasmidIncs[inc].append(metadata[key].ID)
    #plasmidIncs = sorted(plasmidIncs)
    for n in t.traverse(): #loop through the nodes of a tree
        if (n.is_leaf() and n.name == "Reference"):
            #if its the reference branch, populate the faces with column headers
            index = 0

            if len(sensitivePath)>0: #sensitive metadat @ chris
                for sensitive_data_column in sensitive_meta_data.get_columns():
                    (t&"Reference").add_face(addFace(sensitive_data_column), index, "aligned")
                    index = index + 1

            (t&"Reference").add_face(addFace("SampleID"), index, "aligned")
            index = index + 1
            (t&"Reference").add_face(addFace("New?"), index, "aligned")
            index = index + 1
            for i in range(len(plasmidIncs)): #this loop adds the columns (aka the incs) to the reference node
                (t&"Reference").add_face(addFace(list(plasmidIncs.keys())[i]), i + index, "aligned")
            index = index + len(plasmidIncs)
            (t&"Reference").add_face(addFace("MLSTScheme"), index, "aligned")
            index = index + 1
            (t&"Reference").add_face(addFace("Sequence Type"), index, "aligned")
            index = index + 1 
            (t&"Reference").add_face(addFace("Carbapenamases"), index, "aligned")
            index = index + 1
            (t&"Reference").add_face(addFace("Plasmid Best Match"), index, "aligned")
            index = index + 1
            (t&"Reference").add_face(addFace("Best Match Identity"), index, "aligned")
            index = index + 1
            for i in range(len(distanceDict[list(distanceDict.keys())[0]])): #this loop adds the distance matrix
                (t&"Reference").add_face(addFace(distanceDict[list(distanceDict.keys())[0]][i]), index + i, "aligned")
            index = index + len(distanceDict[list(distanceDict.keys())[0]])
        elif (n.is_leaf() and not n.name == "Reference"): 
            #not reference branches, populate with metadata
            index = 0

            if len(sensitivePath)>0: #sensitive metadata @ chris
                # pushing in sensitive data
                for sensitive_data_column in sensitive_meta_data.get_columns():
                    # tree uses bcids like BC18A021A_S12
                    # while sens meta-data uses BC18A021A
                    # trim the "_S.*" if present
                    bcid = str(mData.ID)
                    if bcid.find( "_S" ) != -1:
                        bcid = bcid[ 0:bcid.find( "_S" ) ]
                    sens_col_val = sensitive_meta_data.get_value(bcid=bcid, column_name=sensitive_data_column )
                    n.add_face(addFace(sens_col_val), index, "aligned")
                    index = index + 1

            if (n.name.replace(".fa","") in metadata.keys()):
                mData = metadata[n.name.replace(".fa","")]
            else:
                mData = metadata["na"]
            n.add_face(addFace(mData.ID), index, "aligned")
            index = index + 1
            if (mData.new == True): #new column
                face = e.RectFace(30,30,"green","green") # TextFace("Y",fsize=10,tight_text=True)
                face.border.margin = 5
                face.margin_right = 5
                face.margin_left = 5
                face.vt_align = 1
                face.ht_align = 1
                n.add_face(face, index, "aligned")
            index = index + 1
            for incs in plasmidIncs: #this loop adds presence/absence to the sample nodes
                if (n.name.replace(".fa","") in plasmidIncs[incs]):
                    face = e.RectFace(30,30,"black","black") # TextFace("Y",fsize=10,tight_text=True)
                    face.border.margin = 5
                    face.margin_right = 5
                    face.margin_left = 5
                    face.vt_align = 1
                    face.ht_align = 1
                    n.add_face(face, list(plasmidIncs.keys()).index(incs) + index, "aligned")
            index = index + len(plasmidIncs)
            n.add_face(addFace(mData.MLSTSpecies), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData.SequenceType), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData.CarbapenemResistanceGenes), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData.plasmidBestMatch), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData.plasmididentity), index, "aligned")
            index = index + 1			
            for i in range(len(distanceDict[list(distanceDict.keys())[0]])): #this loop adds distance matrix
                if (n.name in distanceDict): #make sure the column is in the distance matrice
                    n.add_face(addFace(list(distanceDict[n.name])[i]), index + i, "aligned")
                
    t.render(outputFile, w=5000,units="mm", tree_style=ts) #save it as a png, pdf, svg or an phyloxml

    #endregion
#endregion


if __name__ == '__main__':
    start = time.time()#time the analysis
    main()
    end = time.time()
    print("Finished!\nThe analysis used: " + str(end-start) + " seconds")

#!/usr/bin/env python

# This script is a wrapper for CPO Pipeline Phase 3: tree rendering.
# It uses snippy for core genome SNV calling and alignment, clustalw to generate a NJ tree and ete3 to render the dendrogram

import subprocess
import argparse
import csv
import os
import datetime
import errno
import logging
import structlog
import sys
import time
import urllib.request
import uuid
import gzip
import ete3 as e
import drmaa
from pprint import pprint

from pkg_resources import resource_filename

import cpo_pipeline
from cpo_pipeline.pipeline import prepare_job, run_jobs
from cpo_pipeline.tree.parsers import result_parsers


def read(path): #read in a text file to a list
    return [line.rstrip('\n') for line in open(path)]


def addFace(name): #function to add a facet to a tree
    #if its the reference branch, populate the faces with column headers
    face = e.faces.TextFace(name,fsize=10,tight_text=True)
    face.border.margin = 5
    face.margin_right = 5
    face.margin_left = 5
    return face


def main(args, logger=None):
    """
    main entrypoint
    Args:
        args():
    Returns:
        (void)
    """
    
    analysis_id = uuid.uuid4()
    
    curDir = os.getcwd()
    output_dir = args.outdir
    # metadata_file = args.metadata_file
    reference = os.path.abspath(args.reference)
    
    # sensitivePath = str(options.sensitivePath).lstrip().rstrip()
    # sensitiveCols = str(options.sensitiveCols).lstrip().rstrip()
    # outputFile = str(options.outputFile).lstrip().rstrip()
    # bcidCol = str( str(options.bcidCol).lstrip().rstrip() ) 
    # naValue = str( str(options.naValue).lstrip().rstrip() ) 
    
    # metadata = result_parsers.parse_workflow_results(metadata_file)
    # distance = read(distancePath)
    # treeFile = "".join(read(treePath))

    if not logger:
        logging.basicConfig(
            format="%(message)s",
            stream=sys.stdout,
            level=logging.DEBUG,
        )

        structlog.configure_once(
            processors=[
                structlog.stdlib.add_log_level,
                structlog.processors.JSONRenderer()
            ],
            logger_factory=structlog.stdlib.LoggerFactory(),
            wrapper_class=structlog.stdlib.BoundLogger,
            context_class=structlog.threadlocal.wrap_dict(dict),
        )
        logger = structlog.get_logger(
            analysis_id=str(uuid.uuid4()),
            pipeline_version=cpo_pipeline.__version__,
        )

    inputs = []
    with open(args.input_file) as input_file:
        fieldnames = ['sample_id']
        reader = csv.DictReader(
            (row for row in input_file if not row.startswith('#')),
            delimiter='\t',
            fieldnames=fieldnames
        )
        for row in reader:
            inputs.append(row)


    os.environ['QT_QPA_PLATFORM']='offscreen'

    paths = {
        'logs': os.path.abspath(
            os.path.join(
                output_dir,
                'logs',
            )
        ),
        'snippy_output': os.path.abspath(
            os.path.join(
                output_dir,
                "snippy"
            )
        ),
    }

    for output_subdir in paths.values():
        try:
            os.makedirs(output_subdir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    job_script_path = resource_filename('data', 'job_scripts')

    contigs_paths = []
    for sample_id in [input["sample_id"] for input in inputs]:
        contigs = os.path.abspath(
            os.path.join(
                args.result_dir,
                sample_id,
                "assembly",
                "contigs.fa"
            )
        )
        contigs_paths.append(contigs)

    snippy_dirs = [
        os.path.join(
            paths['snippy_output'],
            os.path.basename(
                os.path.dirname(os.path.dirname(contigs)) + "_" + \
                os.path.splitext(os.path.basename(reference))[0]
            )
        )
        for contigs in contigs_paths
    ]
    
    snippy_jobs = [
        {
            'job_name': 'snippy_ctgs',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8 -shell y',
            'remote_command': os.path.join(job_script_path, 'snippy_ctgs.sh'),
            'args': [
                "--ref", reference,
                "--ctgs", contigs,
                "--outdir", snippy_dir
            ]
        }
        for snippy_dir in snippy_dirs
    ]
    
    run_jobs(snippy_jobs)

    
    snippy_core_jobs = [
        {
            'job_name': 'snippy-core',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8 -shell y',
            'remote_command': os.path.join(job_script_path, 'snippy-core.sh'),
            'args': [
                "--ref", reference,
                "--outdir", paths["snippy_output"],
            ] + snippy_dirs
        }
    ]

    run_jobs(snippy_core_jobs)
    
    
    snp_dists_jobs = [
        {
            'job_name': 'snp-dists',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'snp-dists.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.aln"),
            ]
        }
    ]

    run_jobs(snp_dists_jobs)

    iqtree_jobs = [
        {
            'job_name': 'iqtree',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'iqtree.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.full.aln"),
                "--model", "GTR+G4",
            ]
        }
    ]

    run_jobs(iqtree_jobs)


    clonalframeml_jobs = [
        {
            'job_name': 'clonalframeml',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'clonalframeml.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.full.aln"),
                "--treefile", os.path.join(paths["snippy_output"], "core.full.aln.treefile"),
                "--output_file", os.path.join(paths["snippy_output"], "core.full.aln.clonalframeml"),
            ]
        }
    ]

    run_jobs(clonalframeml_jobs)

    maskrc_svg_jobs = [
        {
            'job_name': 'maskrc-svg',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'maskrc-svg.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.full.aln"),
                "--svg", os.path.join(paths["snippy_output"], "core.full.maskrc.svg"),
                "--clonalframeml", os.path.join(paths["snippy_output"], "core.full.aln.clonalframeml"),
                "--output_file", os.path.join(paths["snippy_output"], "core.full.maskrc.aln"),
            ]
        }
    ]
    
    run_jobs(maskrc_svg_jobs)

    snp_sites_jobs = [
        {
            'job_name': 'snp-sites',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'snp-sites.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.full.maskrc.aln"),
                "--output_file", os.path.join(paths["snippy_output"], "core.full.maskrc.snp.aln"),
            ]
        }
    ]

    run_jobs(snp_sites_jobs)

    iqtree_jobs = [
        {
            'job_name': 'iqtree',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'iqtree.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.full.maskrc.aln"),
                "--model", "GTR+G+ASC",
            ]
        }
    ]

    run_jobs(iqtree_jobs)
    
    snp_dists_jobs = [
        {
            'job_name': 'snp-sites',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'snp-dists.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.aln"),
                "--output_file", os.path.join(paths["snippy_output"], "core.matrix.tab"),
            ]
        },
        {
            'job_name': 'snp-sites',
            'output_path': paths['logs'],
            'error_path': paths['logs'],
            'native_specification': '-pe smp 8',
            'remote_command': os.path.join(job_script_path, 'snp-dists.sh'),
            'args': [
                "--alignment", os.path.join(paths["snippy_output"], "core.full.maskrc.snp.aln"),
                "--output_file", os.path.join(paths["snippy_output"], "core.full.maskrc.snp.matrix.tab"),
            ]
        }
    ]

    run_jobs(snp_dists_jobs)
    
    exit(0)
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
        for plasmid in metadata[key]['plasmids']:
            for inc in plasmid['PlasmidRepType'].split(","):
                if (inc.lower().find("inc") > -1):
                    if not (inc in plasmidIncs):
                        plasmidIncs[inc] = [metadata[key]['ID']]
                    else:
                        if metadata[key]['ID'] not in plasmidIncs[inc]:
                            plasmidIncs[inc].append(metadata[key]['ID'])
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

            if (n.name.replace(".fa","") in metadata.keys()):
                mData = metadata[n.name.replace(".fa","")]
            else:
                mData = metadata["na"]
            n.add_face(addFace(mData.ID), index, "aligned")
            index = index + 1
            if (mData['new']): #new column
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
            n.add_face(addFace(mData['MLSTSpecies']), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData['SequenceType']), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData['CarbapenemResistanceGenes']), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData['plasmidBestMatch']), index, "aligned")
            index = index + 1
            n.add_face(addFace(mData['plasmididentity']), index, "aligned")
            index = index + 1			
            for i in range(len(distanceDict[list(distanceDict.keys())[0]])): #this loop adds distance matrix
                if (n.name in distanceDict): #make sure the column is in the distance matrice
                    n.add_face(addFace(list(distanceDict[n.name])[i]), index + i, "aligned")
                
    t.render(outputFile, w=5000,units="mm", tree_style=ts) #save it as a png, pdf, svg or an phyloxml


if __name__ == '__main__':
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name, description='')
    parser.add_argument("-i", "--input", dest="input_file",
                        help="identifier of the isolate", required=True)
    parser.add_argument("-r", "--reference", dest="reference",
                        help="Reference file (fasta or gbk)", required=True)
    parser.add_argument("-m", "--metadata", dest="metadata",
                        help="Metadata file", required=False)
    parser.add_argument("-o", "--outdir", dest="output", default='./',
                        help="absolute path to output folder")
    parser.add_argument("-d", "--resultdir", dest="result_dir",
                        help="Path to cpo-pipeline output folder")
    parser.add_argument('-c', '--config', dest='config_file',
                        default=resource_filename('data', 'config.ini'),
                        help='Config File', required=False)
    args = parser.parse_args()
    main(args)

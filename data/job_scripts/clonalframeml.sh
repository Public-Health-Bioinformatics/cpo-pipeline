#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N clonalframeml
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -t|--treefile TREEFILE -a|--alignment ALIGNMENT -o|--output_file OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

treefile=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -t|--treefile)
    # tree file
    treefile="$2"
    shift # past argument
    shift # past value
    ;;
    -a|--alignment)
    # alignment 
    alignment="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_file)
    # output file
    output_file="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate clonalframeml-1.11

ClonalFrameML \
    "${treefile}" \
    "${alignment}" \
    "${output_file}" \

source deactivate


#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N iqtree
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -a|--alignment aln.fa -m|--model MODEL\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

alignment=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -a|--alignment)
    # alignment file
    alignment="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--model)
    # model
    model="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate iqtree-1.6.7

iqtree \
    -s "${alignment}" \
    -redo \
    -nt AUTO \
    -st DNA \
    -m "${model}" \
    -bb 1000 \
    -alrt 1000 \

source deactivate


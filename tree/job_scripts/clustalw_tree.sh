#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N clustalw
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=2G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -a|--alignment ALIGNMENT \n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

alignment=""
tree_type="nj"

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -a|--alignment)
    # multiple-sequence alignment
    alignment="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate clustalw-2.1

clustalw2 \
    -tree \
    -outputtree="${tree_type}" \
    -infile="${alignment}"

source deactivate

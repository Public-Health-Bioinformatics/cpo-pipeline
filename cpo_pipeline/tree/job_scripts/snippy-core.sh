#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N snippy-core
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=2G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -r|--reference snippy_output_dir1 snippy_output_dir2 ...\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

reference=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -r|--reference)
    # reference genome fasta
    reference="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate snippy-4.2.3

snippy-core \
    --ref "${reference}" \
    $@

source deactivate

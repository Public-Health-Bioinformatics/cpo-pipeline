#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N samtools_faidx
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -f|--fasta\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

fasta=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -f|--fasta)
    # fasta file
    fasta="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate samtools-1.9

samtools faidx "${fasta}"

source deactivate


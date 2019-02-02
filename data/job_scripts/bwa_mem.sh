#!/bin/bash

#$ -V             # Pass environment variables to the job
#$ -N bwa_mem
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
    -1|--R1)
    # input_R1.fastq.gz file
    input_r1_fastq="$2"
    shift # past argument
    shift # past value
    ;;
    -2|--R2)
    # input_R2.fastq.gz file
    input_r2_fastq="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--reference)
    # reference fasta file
    reference="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    # output sam file
    output="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate bwa-0.7.17

bwa mem \
    -t 8 \
    -a \
    -Y \
    -M \
    "${reference}" \
    "${input_r1_fastq}" \
    "${input_r2_fastq}" \
    -o "${output}"

source deactivate


#!/usr/bin/env bash

#$ -V             # Pass environment variables to the job
#$ -N mash_screen
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -1|--R1 INPUT_R1_FASTQ -2|--R2 INPUT_R2_FASTQ -q|--queries QUERIES_MSH -i|--min-identity MIN_IDENTITY -o|--output_file OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

input_r1_fastq=""
input_r2_fastq=""
queries=""
output_file=""
min_identity=0.996

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
    -q|--queries)
    # mash sketch file <queries>.msh
    queries="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--min-identity)
    # Minimum identity to report
    min_identity="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_file)
    # Output file
    output_file="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

mkdir -p $(dirname "${output_file}")

source activate mash-2.0

mash screen \
     -p 8 \
     -w \
     -i "${min_identity}" \
     "${queries}" \
     "${input_r1_fastq}" \
     "${input_r2_fastq}" \
     > "${output_file}"

source deactivate

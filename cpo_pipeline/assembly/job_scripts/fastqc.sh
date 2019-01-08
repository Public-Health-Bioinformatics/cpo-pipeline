#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N fastqc
#$ -cwd           # Use the current working dir
#$ -pe smp 4      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -1|--R1 INPUT_R1_FASTQ -2|--R2 INPUT_R2_FASTQ -o|--output_dir OUTPUT_DIR\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

input_fastq=""
output_dir=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -1|--R1)
    # input.fastq.gz file
    input_r1_fastq="$2"
    shift # past argument
    shift # past value
    ;;
    -2|--R2)
    # input.fastq.gz file
    input_r2_fastq="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_dir)
    # Output directory
    output_dir="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

mkdir -p "${output_dir}"

source activate fastqc-0.11.7

fastqc -t 4 -o "${output_dir}" --extract "${input_r1_fastq}" "${input_r2_fastq}"

source deactivate

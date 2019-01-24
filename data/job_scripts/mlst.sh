#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N mlst
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -l|--label LABEL -i|--input INPUT_CONTIGS_FASTA -o|--output_file OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

input=""
output_file=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -i|--input)
    # input_R1.fastq.gz file
    assembly="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--label)
    # Sample ID
    label="$2"
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

source activate mlst-2.15.1

mlst --label "${label}" "${assembly}" > "${output_file}"

source deactivate

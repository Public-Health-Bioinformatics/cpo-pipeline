#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N seqtk
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -1|--R1 INPUT_R1_FASTQ -2|--R2 INPUT_R2_FASTQ -o|--output_file OUTPUT_FILE\n\
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
output_file=""

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
    -o|--output_file)
    # Output file
    output_file="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

mkdir -p $(dirname "${output_file}")

source activate seqtk-1.3

R1bp=`seqtk fqchk $input_r1_fastq | head -3 | tail -1 | cut -d$'\t' -f 2`
R2bp=`seqtk fqchk $input_r2_fastq | head -3 | tail -1 | cut -d$'\t' -f 2`
totalbp=$((R1bp + R2bp))

echo "${totalbp}" > "${output_file}"

source deactivate


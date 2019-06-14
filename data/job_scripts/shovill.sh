#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N shovill
#$ -cwd           # Use the current working dir
#$ -pe smp 16     # Parallel Environment (how many cores)
#$ -l h_vmem=2G   # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -1|--R1 INPUT_R1_FASTQ -2|--R2 INPUT_R2_FASTQ -c|--mincov MIN_CONTIG_COVERAGE -l|--minlen MIN_CONTIG_LENGTH -o|--output_dir OUTPUT_DIR\n\
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
mincov="2"
minlen="0"
output_dir=""

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
    -c|--mincov)
    # minimum contig coverage
    mincov="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--minlen)
    # minimum contig length
    minlen="$2"
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

source activate shovill-1.0.1

shovill \
    --mincov "${mincov}" \
    --minlen "${minlen}" \
    --force \
    --R1 "${input_r1_fastq}" \
    --R2 "${input_r2_fastq}" \
    --cpus 16 \
    --ram 32 \
    --tmpdir shovill_tmp \
    --outdir "${output_dir}"

source deactivate

rm -rf "${output_dir}"/shovill_tmp

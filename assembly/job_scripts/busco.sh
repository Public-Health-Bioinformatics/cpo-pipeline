#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N busco
#$ -cwd           # Use the current working dir
#$ -pe smp 8     # Parallel Environment (how many cores)
#$ -l h_vmem=4G   # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input INPUT_CONTIGS_FASTA -l|--database BUSCO_DB -o|--output_dir OUTPUT_DIR\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

contigs=""
database=""
output_dir=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -i|--input)
    # input contigs.fa file
    contigs="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--database)
    # busco db
    database="$2"
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

contigs_abspath=$(readlink -f "${contigs}")

cd "${output_dir}"

source activate busco-3.0.2

run_busco \
    --in "${contigs_abspath}" \
    --lineage_path "${database}" \
    --mode genome \
    --out busco \
    --blast_single_core \
    --cpu 8 \
    --force

source deactivate

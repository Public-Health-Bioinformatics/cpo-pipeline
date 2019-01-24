#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N rgi
#$ -cwd           # Use the current working dir
#$ -pe smp 12     # Parallel Environment (how many cores)
#$ -l h_vmem=2G   # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input INPUT_CONTIGS_FASTA -c|--card_json CARD_JSON_FILE -o|--output_dir OUTPUT_DIR\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

assembly=""
card_json=""
output_dir=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -i|--input)
    # input contigs.fa file
    assembly="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--card_json)
    # location of database folders
    card_json="$2"
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

echo "${card_json}"

assembly_realpath=$( realpath "${assembly}" )

mkdir -p "${output_dir}"

echo "going to cd into " "${output_dir}" "..."

cd "${output_dir}"

pwd

source activate rgi-4.0.3

rgi load \
    -i "${card_json}" \
    --local

rgi main \
    -i "${assembly_realpath}" \
    -o rgi \
    -t contig \
    -a BLAST \
    -n 12 \
    --local \
    --clean

rm -rf localDB

source deactivate

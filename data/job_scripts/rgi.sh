#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N rgi
#$ -cwd           # Use the current working dir
#$ -pe smp 12     # Parallel Environment (how many cores)
#$ -l h_vmem=2G   # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input INPUT_CONTIGS_FASTA -c|--card_json CARD_JSON_FILE -o|--output_file OUTPUT_FILE\n\
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
output_file=""

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
    -o|--output_file)
    # Output file
    output_file="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

echo "${card_json}"

mkdir -p $(dirname "${output_file}")

echo "going to cd into " $(dirname "${output_file}") "..."

cd $(dirname "${output_file}")

pwd

source activate rgi-4.0.3

rgi load \
    -i "${card_json}" \
    --local

rgi main \
    -i "${assembly}" \
    -o $(basename "${output_file}") \
    -t contig \
    -a BLAST \
    -n 12 \
    --local \
    --clean

rm -rf localDB

source deactivate

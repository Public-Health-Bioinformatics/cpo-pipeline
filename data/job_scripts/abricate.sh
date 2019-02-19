#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N abricate
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=2G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] [-d|--datadir ABRICATE_DATA_DIR] -i|--input INPUT_CONTIGS_FASTA -D|--database ABRICATE_DATABASE_NAME -o|--output_file OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

assembly=""
datadir=""
db=""
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
    -d|--datadir)
    # location of database folders
    datadir="$2"
    shift # past argument
    shift # past value
    ;;
    -D|--database)
    # database to use
    db="$2"
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

source activate abricate-0.8.7

abricate \
    $( if [ "${datadir}" != "" ]; then echo "--datadir" "${datadir}"; fi ) \
    --db "${db}" \
    "${assembly}" > \
    "${output_file}"

source deactivate

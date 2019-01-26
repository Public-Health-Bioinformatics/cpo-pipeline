#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N quast
#$ -cwd           # Use the current working dir
#$ -pe smp 8     # Parallel Environment (how many cores)
#$ -l h_vmem=4G   # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] [-r|--reference REFERENCE_FASTA] -i|--input INPUT_CONTIGS_FASTA -o|--outdir OUTPUT_DIR\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

contigs=""
outdir=""
reference=""

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
    -r|--reference)
    # reference fasta file
    reference="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outdir)
    # Output directory
    outdir="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

mkdir -p "${outdir}"

source activate quast-5.0.2

quast \
    --output-dir "${outdir}" \
    --threads 8 \
    --fast \
    --silent \
    --conserved-genes-finding \
    $( if [ "${reference}" != "" ]; then echo "-r" "${reference}"; fi ) \
    "${contigs}" \

source deactivate

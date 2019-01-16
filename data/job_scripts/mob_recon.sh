#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N mob_recon
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=2G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input INPUT_CONTIGS_FASTA -o|--output_dir OUTPUT_DIR\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

assembly=""
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
    -o|--output_dir)
    # Output file
    output_dir="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

mkdir -p "${output_dir}"

source activate mob_suite-1.4.5

#run mob typer
mob_recon \
    --infile "${assembly}" \
    --outdir "${output_dir}" \
    --run_typer

source deactivate

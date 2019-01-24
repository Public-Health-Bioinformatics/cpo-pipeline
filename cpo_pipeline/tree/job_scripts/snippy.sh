#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N snippy
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=2G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -r|--reference REFERENCE_GENOME -c|--contigs CONTIGS_FASTA -o|--output_dir OUTPUT_DIR\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

reference=""
contigs=""=""
output_dir=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -r|--reference)
    # reference genome
    reference="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--ctgs)
    # contigs.fa file
    contigs="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_dir)
    # output directory
    output_dir="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate snippy-4.2.3

snippy \
    --cpus 8 \
    --reference "${reference}" \
    --ctgs "${contigs}" \
    --outdir "${output_dir}"  

source deactivate

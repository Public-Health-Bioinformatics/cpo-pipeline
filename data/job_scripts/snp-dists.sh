#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N snp-dists
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -a|--alignment aln.fa  -o|--output_file OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

alignment=""
output_file=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -a|--alignment)
    # alignment fasta file
    alignment="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_file)
    # Output directory
    output_file="$2"
    shift # past argument
    shift # past value
    ;;

  esac
done

mkdir -p `dirname "${output}"`

source activate snp-dists-0.6.3

snp-dists \
    -a "${alignment}" \
    > "${output_file}"

source deactivate


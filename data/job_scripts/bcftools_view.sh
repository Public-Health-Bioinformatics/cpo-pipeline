#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N bcftools_view
#$ -cwd           # Use the current working dir
#$ -pe smp 2      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input VCF -o|--output VCF \n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

input=""
output=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -i|--input)
    # input sam (or bam) file
    input="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    # output file
    output="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate bcftools-1.9

bcftools view \
	 -Ov \
	 -v snps,mnps \
	 -o "${output}" \
	 "${input}"

source deactivate


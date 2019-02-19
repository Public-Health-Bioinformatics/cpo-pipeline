#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N freebayes
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input SAM_BAM -r|--reference REFERENCE_FASTA -o|--output VCF \n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

threads=8

input=""
output=""
reference=""

ploidy=1
min_base_quality=20
min_mapping_quality=60
min_coverage=10
min_alternate_fraction=0.8
min_repeat_entropy=1.0

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
    -r|--reference)
    # reference fasta
    reference="$2"
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

reference_length=$( cut -f 2 "${reference}".fai )
region_length=$(( $reference_length / $threads ))

source activate freebayes-1.2

freebayes-parallel \
    <(fasta_generate_regions.py "${reference}" "${region_length}") \
    "${threads}" \
    --ploidy 1 \
    --min-base-quality "${min_base_quality}" \
    --min-mapping-quality "${min_mapping_quality}" \
    --min-coverage "${min_coverage}" \
    --min-alternate-fraction "${min_alternate_fraction}" \
    --min-repeat-entropy "${min_repeat_entropy}" \
    -f "${reference}" \
    "${input}" \
    > "${output}"

source deactivate


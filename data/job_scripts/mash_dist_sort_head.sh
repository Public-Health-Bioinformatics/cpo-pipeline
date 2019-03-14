#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N mash_dist
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] [-m MINIMUM_KMERS] -1|--R1 INPUT_R1_FASTQ -2|--R2 INPUT_R2_FASTQ -q|--queries QUERIES_MSH -o|--output_file OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

input_r1_fastq=""
input_r2_fastq=""
queries=""
minimum_kmers=3
output_file=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -1|--R1)
    # input_R1.fastq.gz file
    input_r1_fastq="$2"
    shift # past argument
    shift # past value
    ;;
    -2|--R2)
    # input_R2.fastq.gz file
    input_r2_fastq="$2"
    shift # past argument
    shift # past value
    ;;
    -q|--queries)
    # mash sketch file <queries>.msh
    queries="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--minimum-kmers)
    # ominimum kmers
    minimum_kmers="$2"
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

output_dir=$( dirname "${output_file}" )
mkdir -p "${output_dir}"

source activate mash-2.0

cat "${input_r1_fastq}" "${input_r2_fastq}" | mash sketch -m "${minimum_kmers}" -r -o "${output_dir}"/reads -
mash dist -p 8 "${queries}" "${output_dir}"/reads.msh | sort -gk2,2 | head > "${output_file}"
rm "${output_dir}"/reads.msh

source deactivate

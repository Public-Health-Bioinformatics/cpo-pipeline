#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N samtools_sort
#$ -cwd           # Use the current working dir
#$ -pe smp 4      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -i|--input SAM_BAM -o|--output SAM_BAM -n|--name-order\n\
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
name_order=false

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
    -n|--name-order)
    # sort by readname
    name_order=true
    shift # past argument
    ;;
    -o|--output)
    # output sam (or bam) file
    output="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate samtools-1.9

samtools sort \
	 --threads 3 \
	 $( if [ "$name_order" = true ]; then printf "%s" "-n"; fi ) \
	 -l 0 \
	 "${input}" \
	 -o "${output}"

source deactivate


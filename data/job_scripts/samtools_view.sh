#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N samtools_view
#$ -cwd           # Use the current working dir
#$ -pe smp 4      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] [-F|--flags FLAGS] -i|--input SAM_BAM|\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

input=""
flags=0

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
    -F|--flags)
    # only include reads with none of the FLAGS in this integer present
    flags="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate samtools-1.9

samtools view \
	 -@ 3 \
	 -h \
	 -F "${flags}" \
	 "${input}"
	 

source deactivate


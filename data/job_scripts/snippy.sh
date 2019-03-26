#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N snippy
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -r|--ref REFERENCE_FASTA -1|--R1 INPUT_R1_FASTQ -2|--R2 INPUT_R2_FASTQ  -o|--outdir OUTPUT_FILE\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

reference=""
input_r1_fastq=""
input_r2_fastq=""
outdir=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -r|--ref)
    # reference fasta file
    reference="$2"
    shift # past argument
    shift # past value
    ;;
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
    -o|--outdir)
    # Output directory
    outdir="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

mkdir -p `dirname "${outdir}"`

source activate snippy-4.3.6

snippy \
    --ref "${reference}" \
    --R1 "${input_r1_fastq}" \
    --R2 "${input_r2_fastq}" \
    --cpus 8 \
    --ram 8 \
    --outdir "${outdir}" \
    --cleanup

source deactivate


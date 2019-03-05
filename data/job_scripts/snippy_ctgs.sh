#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N snippy_contigs
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -r|--ref REFERENCE_FASTA -c|--ctgs -o|--outdir OUTPUT_FILE\n\
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
    -c|--ctgs)
    # contigs fasta file
    ctgs="$2"
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
    --ctgs "${ctgs}" \
    --cpus 8 \
    --ram 8 \
    --outdir "${outdir}"

source deactivate


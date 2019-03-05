#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N snippy-core
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -r|--ref REFERENCE_FASTA  -o|--outdir OUTPUT_FILE SNIPPY_DIR [SNIPPY_DIR ...]\n\
\n\
optional arguments:\n\
  -h, --help \t\t\t Show this help message and exit" 

if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
  echo -e "Usage: ${USAGE}"
  exit 0
fi

reference=""
outdir=""
POSITIONAL=()

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
    -o|--outdir)
    # Output directory
    outdir="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

mkdir -p "${outdir}"

cd "${outdir}"

source activate snippy-4.3.6

snippy-core \
    --ref "${reference}" \
    "$@"

source deactivate


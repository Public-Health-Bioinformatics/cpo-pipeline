#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N maskrc-svg
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=1G   # Memory (RAM) allocation *per core*

USAGE="qsub $( basename "$BASH_SOURCE" ) [-h] -a|--alignment ALIGNMENT -s|--svg SVG_OUTPUT -o|--output_file OUTPUT_FILE -c|--clonalframeml CLONALFRAMEML_PREFIX\n\
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
clonalframeml=""
svg=""

while [[ $# -gt 0 ]]
do
  key="$1"
  
  case $key in
    -a|--alignment)
    # alignment 
    alignment="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--clonalframeml)
    # prefix used for clonalframeml output
    clonalframeml="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--svg)
    # svg image output file
    svg="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_file)
    # output file
    output_file="$2"
    shift # past argument
    shift # past value
    ;;
  esac
done

source activate maskrc-svg-0.5

maskrc-svg.py  \
    --aln "${alignment}" \
    --out "${output_file}" \
    --svg "${svg}" \
    "${clonalframeml}" \

source deactivate


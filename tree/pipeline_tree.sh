#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

#input parameters: 

source activate snippy-4.2.3

threads=8
ram=80
outputDir="$1"
reference="$2"
contigs="$3"

echo "parameters: "
echo "outputDir: $outputDir"
echo "reference: $reference"
echo "contigs: $contigs"
echo "threads: $threads"
echo "ram: $ram"

tree_dir="$outputDir"/tree
mkdir -p "${tree_dir}"

contigArray=()

#find snps within each genome
IFS=',' read -ra ADDR <<< "$contigs" #hax to read in a csv
for i in "${ADDR[@]}"; do
	refG=`basename $i`
	echo "${refG}"
	snippy --cpus "${threads}" --outdir "${tree_dir}"/"${refG}" --reference "${reference}" --ctgs $i
	contigArray+=("${refG}")
done

cd "${tree_dir}"

#create an alignment from the snps
snippy-core --prefix core *

#create a distance matrix from alignment
snp-dists core.full.aln > distance.tab

#make a nj tree from the matrix
clustalw2 -tree -infile=core.full.aln -outputtree=nj

source deactivate
echo "done creating Neighbour Joining tree"





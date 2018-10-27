#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

#input parameters: 

threads=8
ram=80
output_dir="$1"
reference="$2"
contig_dir="$3"

echo "parameters: "
echo "output_dir: $output_dir"
echo "reference: $reference"
echo "contig_dir: $contig_dir"
echo "threads: $threads"
echo "ram: $ram"

tree_dir="${output_dir}"/tree
mkdir -p "${tree_dir}"

contigArray=()

source activate snippy-4.2.3

#find snps within each genome

for contig in "${contig_dir}"/*.fa; do
	refG=`basename "${contig}"`
	echo "${refG}"
	snippy --cpus "${threads}" --outdir "${tree_dir}"/"${refG}" --reference "${reference}" --ctgs $contig
	contigArray+=("${refG}")
done

cd "${tree_dir}"

#create an alignment from the snps
snippy-core --ref "${reference}" --prefix core *

source deactivate

source activate snp-dists-0.6.2

#create a distance matrix from alignment
snp-dists core.full.aln > distance.tab

source deactivate

source activate clustalw-2.1

clustalw2 -tree -infile=core.full.aln -outputtree=nj

source deactivate

echo "done creating Neighbour Joining tree"





#!/bin/bash

fastq_dir="/projects/carbapenamase_producing_organisms_surveillance/fastq_symlinks"
sample_list=$1
declare -a reads=("R1" "R2")

for sample_id in `cat "${sample_list}"`; do
    for read in "${reads[@]}"; do
	if [ -f "${fastq_dir}"/"${sample_id}"_S*_"${read}"_*.fastq.gz ]; then
	    ln -s "${fastq_dir}"/"${sample_id}"_S*_"${read}"_*.fastq.gz "${sample_id}"_"${read}".fastq.gz
	else
	    (>&2 echo "${sample_id}"_"${read}" not found.)
	fi
    done
done

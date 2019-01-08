#!/bin/bash

#parameters:
# /home/dfornika/code/CPO_Prediction/pipeline.qsub
# /home/dfornika/code/CPO_Prediction/pipeline.py
# samples.txt
# fastq_input_dir
# output_dir
# species_name ("Escherichia coli")

pipeline_qsub_script=$1
pipeline_script=$2
samples_list=$3
input_dir=$4
output_dir=$5
species_name=$6

#submitting multiple jobs
for sample_id in `cat "${samples_list}"`; do

    reads_1="${input_dir}"/"${sample_id}"_R1.fastq.gz
    reads_2="${input_dir}"/"${sample_id}"_R2.fastq.gz
    echo "${sample_id}"
    echo "${reads_1}"
    echo "${reads_2}"
    
    qsub  "${pipeline_qsub_script}" \
	  "${pipeline_script}" \
	  "${sample_id}" \
	  "${reads_1}" \
	  "${reads_2}" \
	  "${output_dir}" \
	  "${species_name}";
done;

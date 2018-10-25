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
module_1_output_dir=$4
module_2_output_dir=$5

#submitting multiple jobs
for sample_id in `cat "${samples_list}"`; do

    assembly="${module_1_output_dir}"/contigs/"${sample_id}".fa
    echo "${sample_id}"
    
    qsub  "${pipeline_qsub_script}" \
	  "${pipeline_script}" \
	  "${sample_id}" \
	  "${assembly}" \
	  "${module_2_output_dir}"
done;

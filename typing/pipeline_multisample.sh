#!/bin/bash

#parameters:
# /home/dfornika/code/CPO_Prediction/pipeline.qsub
# /home/dfornika/code/CPO_Prediction/pipeline.py
# samples.txt
# fastq_input_dir
# output_dir
# species_name ("Escherichia coli")

pipeline_qsub_script=$1
samples_list=$2
module_1_output_dir=$3
module_2_output_dir=$4

#submitting multiple jobs
for sample_id in `cat "${samples_list}"`; do

    assembly="${module_1_output_dir}"/contigs/"${sample_id}".fa
    echo "${sample_id}"
    
    qsub  "${pipeline_qsub_script}" \
	  "${sample_id}" \
	  "${assembly}" \
	  "${module_2_output_dir}" \
	  "/data/ref_databases/card/card-2.0.3.json" \
	  "/projects/carbapenemase_producing_organisms_surveillance/databases/abricate" \
	  "carbapenemases"
	  
done;

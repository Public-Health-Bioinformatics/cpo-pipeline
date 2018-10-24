#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

################################################################################
# J-J-J-Jia @ pipeline_assembly_contaminant.sh: filter read k-mers to reference genome(s) then assemble the filtered reads. then do QA on them with metaquast and busco
#
# input parameters: 1=sample_id, 2=forward, 3 = reverse, 4=output, 5=tmpdir (for shovill), 6=reference genomes (csv, no spaces), 7=buscodb, 8=correct reference
#
# Requires: shovill, bbmap, quast and busco (all conda-ble)
#
# pipeline_assembly_contaminant.sh BC11_Kpn05 BC11_Kpn05_R1.fastq.gz BC11_Kpn05_R2.fastq.gz output /tmp/shovill "ref1.fna,ref2.fna,ref3.fna" enterobacteriales_odb9
################################################################################

sample_id="$1"
reads1_file="$2"
reads2_file="$3"
output_dir="$4"
tmp_dir="$5"
reference_genome="$6"
busco_db="$7"
correct_reference="$8"
threads=8
ram=80


echo "parameters: "
echo "sample_id: ${sample_id}"
echo "reads1_file: ${reads1_file}"
echo "reads2_file: ${reads2_file}"
echo "output_dir: ${output_dir}"
echo "tmp_dir: ${tmp_dir}"
echo "reference_genome: ${reference_genome}"
echo "busco_db: ${busco_db}"
echo "correct_refernce: ${correct_reference}"
echo "threads: $threads"
echo "ram: $ram"


assembly_dir="${output_dir}"/assembly
contigs_dir="${output_dir}"/contigs
temp_dir="${tmp_dir}"/"${sample_id}"
split_output_dir="${assembly_dir}"/filteredReads
qc_dir="${output_dir}"/assembly_qc

#step2, shovill assembly

mkdir -p "${assembly_dir}"
mkdir -p "${contigs_dir}"
mkdir -p "${temp_dir}"
mkdir -p "${split_output_dir}"
mkdir -p "${qc_dir}"


echo "step2: assembly"

#mkdir -p "${assembly_dir}/${sample_id}"

source activate bbmap-38.22

#try to seperate out the contaminant reads
cd "${temp_dir}" #bbsplit creates temp files in the current directory. That cannot be changed. cd first to keep it clean
bbsplit.sh in1="${reads1_file}" in2="${reads2_file}" ref="${reference_genome}" \
	   basename="${split_output_dir}"/"${sample_id}"%_#.fq \
	   outu1="${split_output_dir}"/"${sample_id}"_unmap1.fq \
	   outu2="${split_output_dir}"/"${sample_id}"_unmap2.fq \
	   t="$threads" -Xmx"$ram"g

source deactivate

#assembled all the contaminant reads
IFS=',' read -ra ADDR <<< "${reference_genome}" #hax to read in a csv
for i in "${ADDR[@]}"; do
	refG=`basename $i`
	echo "${refG}"
	assembly_output_dir="${assembly_dir}/${sample_id}/${refG}"
	reads1_proper="${split_output_dir}"/"${sample_id}""${refG}"_1.fq
	reads2_proper="${split_output_dir}"/"${sample_id}""${refG}"_2.fq

	source activate shovill-1.0.1
	
	shovill --R1 "${reads1_proper}" --R2 "${reads2_proper}" \
		--cpus "$threads" --ram "$ram" --tmpdir "${temp_dir}" \
		--outdir "${assembly_output_dir}"

	source deactivate

	contig_proper="${contigs_dir}"/"${sample_id}"."${refG}".fa
	
	# Move assemblies to $assembly_output_dir
	if [ -f "${assembly_output_dir}/contigs.fa" ]
	then
		cp "${assembly_output_dir}/contigs.fa" "${contig_proper}"
	else
		echo "!!!Error during assembly"
		exit 818
	fi

	mkdir -p "${qc_dir}"/"${sample_id}"
	cd "${qc_dir}"/"${sample_id}"

	#run metaquast on 1 of the contaminant genomes
	
	source activate quast-4.6.3
	
	metaquast "${contig_proper}" -R "${reference_genome}" --threads "$threads" \
		  -o "${qc_dir}/${sample_id}/${sample_id}.${refG}.quast"
	
	source deactivate
done

cd "${qc_dir}"/"${sample_id}"

# Run busco on all the assembled contaminant genomes

source activate busco-3.0.2

IFS=',' read -ra ADDR <<< "${reference_genome}" #hax to read in a csv
for i in "${ADDR[@]}"; do
	refG=`basename $i`
	echo "${refG}"
	contig_proper="${contigs_dir}"/"${sample_id}"."${refG}".fa
	#have to run busco with 1 threads cuz it throws some error about tblastn crashing when multithreadding
	run_busco -i "${contig_proper}" -o "${sample_id}.${refG}.busco" \
		  -l "${busco_db}" -m genome -c 1 -f
	mv "run_${sample_id}.${refG}.busco" "${sample_id}.${refG}.busco"
done

source deactivate

# Rename the correct assembly and qc files
mv "${contigs_dir}"/"${sample_id}"."${correct_reference}".fa "${contigs_dir}"/"${sample_id}".fa
mv "${qc_dir}"/"${sample_id}"/"${sample_id}"."${correct_reference}".busco "${qc_dir}"/"${sample_id}"/"${sample_id}".busco
mv "${qc_dir}"/"${sample_id}"/"${sample_id}"."${correct_reference}".quast "${qc_dir}"/"${sample_id}"/"${sample_id}".quast

# Cleanup
rm -rf "${qc_dir}"/"${sample_id}"/tmp
rm -rf "${temp_dir}"

echo "done step 2"

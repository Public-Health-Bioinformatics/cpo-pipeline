#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

######################################################################
# J-J-J-Jia @ pipeline_assembly.sh: assemble the reads. Then do QA on them with busco and quast.
# 
# Positional input parameters:
# 1 = id
# 2 = forward
# 3 = reverse
# 4 = output
# 5 = tmpdir for shovill
# 6 = reference genome
# 7 = buscoDB
# 
# Requires: shovill, bbmap, quast and busco
# 
# Example usage:
# pipeline_assembly.sh BC11 BC11-Kpn005_S2_L001_R1_001.fastq.gz BC11-Kpn005_S2_L001_R2_001.fastq.gz output /home/dfornika/tmp/shovill /home/jjjjia/testCases/tests/references/refbc11 /data/ref_databases/busco/enterobacteriales_odb9
######################################################################

sample_id="$1"
reads1_file="$2"
reads2_file="$3"
output_dir="$4"
tmp_dir="$5"
reference_genome="$6"
busco_db="$7"
threads=8
ram=80

echo "parameters: "
echo "sample_id: ${sample_id}"
echo "reads1_file: ${reads1_file}"
echo "reads2_file: ${reads2_file}"
echo "output_dir: ${output_dir}"
echo "tmp_dir: ${tmp_dir}"
echo "reference_genome: ${ref_genome}"
echo "busco_db: ${busco_db}"
echo "threads: $threads"
echo "ram: $ram"

assembly_dir="${output_dir}"/assembly
contigs_dir="${output_dir}"/contigs
qc_dir="${output_dir}"/assembly_qc
temp_dir="${tmp_dir}"/"${sample_id}"

#step2, shovill assembly

mkdir -p "${assembly_dir}"
mkdir -p "${contigs_dir}"
mkdir -p "${qc_dir}"
mkdir -p "${temp_dir}"

echo "step2: assembly"
assembly_output_dir="${assembly_dir}/${sample_id}"

source activate shovill-1.0.1

shovill \
    --mincov 3 \
    --minlen 500 \
    --force \
    --R1 "${reads1_file}" \
    --R2 "${reads2_file}" \
    --cpus "${threads}" \
    --ram "${ram}" \
    --tmpdir "${temp_dir}" \
    --outdir "${assembly_output_dir}"

source deactivate

# Move assemblies to $assembly_output_dir
if [ -f "${assembly_output_dir}/contigs.fa" ]
then
	cp "${assembly_output_dir}/contigs.fa" "${contigs_dir}/${sample_id}.fa"
else
	echo "!!!Error during assembly"
	exit 818
fi

source activate quast-4.6.3

#run quast on assembled genome
mkdir -p "${qc_dir}"/"${sample_id}"
quast \
    "${contigs_dir}/${sample_id}.fa" \
    -R "${reference_genome}" \
    -o "${qc_dir}/${sample_id}/${sample_id}.quast" \
    --threads "${threads}"

source deactivate

source activate busco-3.0.2

cd "${qc_dir}"/"${sample_id}"
run_busco \
    -i "../../../${contigs_dir}/${sample_id}.fa" \
    -o "${sample_id}.busco" \
    -l "${busco_db}" \
    -m genome \
    -c 1 \
    -sp E_coli_K12 \
    -f

mv "run_${sample_id}".busco "${sample_id}".busco

source deactivate

rm -rf "${qc_dir}"/"${sample_id}"/tmp
rm -rf "${temp_dir}"

echo "done step 2"

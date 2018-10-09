#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N CPO_pipeline    # Replace with a more specific job name
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

#####################################################################################################################################################################################################
#J-J-J-Jia @ pipeline_assembly.sh: assemble the reads.	then do QA on them	with busco and quast																									#
#input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome, 7=buscoDB																					#
#Requires: shovill, bbmap, quast and busco (all conda-ble)																																			#
# pipeline_assembly.sh BC11 BC11-Kpn005_S2_L001_R1_001.fastq.gz BC11-Kpn005_S2_L001_R2_001.fastq.gz output /home/dfornika/tmp/shovill /home/jjjjia/testCases/tests/references/refbc11 /data/ref_databases/busco/enterobacteriales_odb9
#####################################################################################################################################################################################################

ID="$1"
R1="$2"
R2="$3"
outputDir="$4"
tmpDir="$5"
refGenome="$6"
buscoDB="$7"
threads=8
ram=80

echo "parameters: "
echo "ID: $ID"
echo "R1: $R1"
echo "R2: $R2"
echo "outputDir: $outputDir"
echo "tmpDir: $tmpDir"
echo "refGenomePath: $refGenome"
echo "buscoDB path: $buscoDB"
echo "threads: $threads"
echo "ram: $ram"

assemblyDir="$outputDir"/assembly
contigsDir="$outputDir"/contigs
qcDir="$outputDir"/assembly_qc
tempDir="$tmpDir"/"$ID"

#step2, shovill assembly

mkdir -p "$assemblyDir"
mkdir -p "$contigsDir"
mkdir -p "$qcDir"
mkdir -p "$tempDir"

echo "step2: assembly"
assemblyOutDir="$assemblyDir/$ID"
#mkdir -p "$assemblyDir/$ID"

source activate shovill-1.0.1

shovill --mincov 3 --minlen 500 --force --R1 "$R1" --R2 "$R2" --cpus "$threads" --ram "$ram" --tmpdir "$tempDir" --outdir "$assemblyOutDir"
#make the contigs dir in cwd

source deactivate

#move all the assemblies to the new phone.
if [ -f "$assemblyOutDir/contigs.fa" ]
then
	cp "$assemblyOutDir/contigs.fa" "$contigsDir/$ID.fa"
else
	echo "!!!Error during assembly"
	exit 818
fi

source activate quast-4.6.3

#run quast on assembled genome
mkdir -p "$qcDir"/"$ID"
quast "$contigsDir/$ID.fa" -R "$refGenome" -o "$qcDir/$ID/$ID.quast" --threads "$threads"

source deactivate

#contamination genomes
#bbsplit.sh in1=BC16-Cfr035_S10_L001_R1_001.fastq.gz in2=BC16-Cfr035_S10_L001_R2_001.fastq.gz ref=ref1,ref2,ref3 basename=o%_#.fq outu1=unmap1.fq outu2=unmap2.fq
#metaquast BC16-Cfr035_S10.fa -R ref1,ref2,ref3 --threads 8 -o result

source activate busco-3.0.2

cd "$qcDir"/"$ID"
run_busco -i "../../../$contigsDir/$ID.fa" -o "$ID.busco" -l $buscoDB -m genome -c 1 -sp E_coli_K12 -f
mv "run_$ID".busco "$ID".busco

source deactivate

rm -rf "$qcDir"/"$ID"/tmp
rm -rf "$tempDir"

echo "done step 2"

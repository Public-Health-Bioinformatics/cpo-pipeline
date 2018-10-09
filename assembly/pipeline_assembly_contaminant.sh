#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

#####################################################################################################################################################################################################
#J-J-J-Jia @ pipeline_assembly_contaminant.sh: filter read k-mers to reference genome(s) then assemble the filtered reads.	then do QA on them	with metaquast and busco							#
#input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome (csv, no spaces), 7=buscodb, 8=correct reference											#		
#Requires: shovill, bbmap, quast and busco (all conda-ble)																																			#
#./pipeline_assembly_contaminant.sh BC11_Kpn05 ~/testCases/seqs/R1/BC11_Kpn05_R1.fastq.gz ~/testCases/seqs/R2/BC11_Kpn05_R2.fastq.gz ~/cpo_assembly ~/cpo_temp "~/ref1.fna,~/ref2.fna,~/ref3.fna"	#
#####################################################################################################################################################################################################
#~/scripts/pipeline_assembly_contaminant.sh BC16 /data/jjjjia/R1/BC16-Cfr035_S10_L001_R1_001.fastq.gz /data/jjjjia/R2/BC16-Cfr035_S10_L001_R2_001.fastq.gz /home/jjjjia/testCases/tests /home/jjjjia/testCases/tests/shovilltemp /home/jjjjia/testCases/tests/references/ref1,/home/jjjjia/testCases/tests/references/ref2,/home/jjjjia/testCases/tests/references/ref3 /home/jjjjia/databases/enterobacteriales_odb9/ 

ID="$1"
R1="$2"
R2="$3"
outputDir="$4"
tmpDir="$5"
refGenome="$6"
buscoDB="$7"
correctRef="$8"
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
echo "correctRefernce: $correctRef"
echo "threads: $threads"
echo "ram: $ram"


assemblyDir="$outputDir"/assembly
contigsDir="$outputDir"/contigs
tempDir="$tmpDir"/"$ID"
splitOutDir="$assemblyDir"/filteredReads
qcDir="$outputDir"/assembly_qc

#step2, shovill assembly

mkdir -p "$assemblyDir"
mkdir -p "$contigsDir"
mkdir -p "$tempDir"
mkdir -p "$splitOutDir"
mkdir -p "$qcDir"


echo "step2: assembly"

#mkdir -p "$assemblyDir/$ID"

source activate bbmap-38.22

#try to seperate out the contaminant reads
cd $tempDir #bbsplit creates temp files in the current directory that cannot be changed. cd first to keep it clean
bbsplit.sh in1="$R1" in2="$R2" ref="$refGenome" basename="$splitOutDir"/"$ID"%_#.fq outu1="$splitOutDir"/"$ID"_unmap1.fq outu2="$splitOutDir"/"$ID"_unmap2.fq t="$threads" -Xmx"$ram"g

source deactivate

#assembled all the contaminant reads
IFS=',' read -ra ADDR <<< "$refGenome" #hax to read in a csv
for i in "${ADDR[@]}"; do
	refG=`basename $i`
	echo $refG
	assemblyOutDir="$assemblyDir/$ID/$refG"
	R1Proper="$splitOutDir"/"$ID""$refG"_1.fq
	R2Proper="$splitOutDir"/"$ID""$refG"_2.fq

	source activate shovill-1.0.1
	
	shovill --R1 "$R1Proper" --R2 "$R2Proper" --cpus "$threads" --ram "$ram" --tmpdir "$tempDir" --outdir "$assemblyOutDir" 

	source deactivate
	
	contigProper="$contigsDir"/"$ID"."$refG".fa
	#move all the assemblies to the new phone.
	if [ -f "$assemblyOutDir/contigs.fa" ]
	then
		cp "$assemblyOutDir/contigs.fa" "$contigProper"
	else
		echo "!!!Error during assembly"
		exit 818
	fi

	mkdir -p "$qcDir"/"$ID"
	cd "$qcDir"/"$ID"

	#run metaquast on 1 of the contaminant genomes
	
	source activate quast-4.6.3
	
	metaquast "$contigProper" -R "$refGenome" --threads "$threads" -o "$qcDir/$ID/$ID.$refG.quast"
	
	source deactivate
done

cd "$qcDir"/"$ID"

#run busco on all the assembled contaminant genomes

source activate busco-3.0.2

IFS=',' read -ra ADDR <<< "$refGenome" #hax to read in a csv
for i in "${ADDR[@]}"; do
	refG=`basename $i`
	echo $refG
	contigProper="$contigsDir"/"$ID"."$refG".fa
	#have to run busco with 1 threads cuz it throws some error about tblastn crashing when multithreadding
	run_busco -i "$contigProper" -o "$ID.$refG.busco" -l "$buscoDB" -m genome -c 1 -f
	mv "run_$ID.$refG.busco" "$ID.$refG.busco"
done

source deactivate

#rename the correct assembly and qc files
mv "$contigsDir"/"$ID"."$correctRef".fa "$contigsDir"/"$ID".fa
mv "$qcDir"/"$ID"/"$ID"."$correctRef".busco "$qcDir"/"$ID"/"$ID".busco 
mv "$qcDir"/"$ID"/"$ID"."$correctRef".quast "$qcDir"/"$ID"/"$ID".quast 

#cleanup
rm -rf cd "$qcDir"/"$ID"/tmp
rm -rf "$tempDir"

echo "done step 2"

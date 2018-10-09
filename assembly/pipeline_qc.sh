#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N CPO_pipeline    # Replace with a more specific job name
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

#####################################################################################################################################################################################################
#J-J-J-Jia @ pipeline_qc.sh: runs seqtk, mash, kraken2 and fastqc for pre-assembly quality checks																									#
#input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=mashgenomerefdb, $6=mashplasmidrefdb, $7=kraken2db, $8=kraken2plasmiddb															#		
#Requires: mash, kraken2, fastqc, seqtk	(all conda-ble)																																				#
# pipeline_qc.sh BC16-Cfr035 BC16-Cfr035_S10_L001_R1_001.fastq.gz BC16-Cfr035_S10_L001_R2_001.fastq.gz output /data/ref_databases/mash/refseq.genomes.k21s1000.msh /data/ref_databases/mash/refseq.plasmid.k21s1000.msh /data/ref_databases/kraken2/2018-09-20_standard /data/ref_databases/kraken2/2018-09-20_plasmid	#
#####################################################################################################################################################################################################

#step 1, mash QC

ID="$1"
R1="$2"
R2="$3"
outputDir="$4"
refDB="$5"
plasmidRefDB="$6"
k2db="$7"
k2plasmid="$8"
threads=8
ram=80

echo "parameters: "
echo "ID: $ID"
echo "R1: $R1"
echo "R2: $R2"
echo "outputDir: $outputDir"
echo "refDB: $refDB"
echo "plasmidDB: $plasmidRefDB"
echo "kraken2db: $k2db"
echo "kraken2plasmiddb: $k2plasmid"
echo "threads: $threads"
echo "ram: $ram"

mashResultDir="$outputDir"/qcResult

mkdir -p "$outputDir"
mkdir -p "$mashResultDir"
mkdir -p "$outputDir"/summary


#step1, mash qc
echo "step1: mash"
qcOutDir="$mashResultDir/$ID"
mkdir -p "$mashResultDir"/"$ID"

cat "$R1" "$R2" > "$qcOutDir"/concatRawReads.fastq
cat "$R1" > "$qcOutDir"/R1.fastq

source activate mash-2.0

#get estimation of genome size (k-mer method)
mash sketch -m 3 "$qcOutDir"/R1.fastq -o "$qcOutDir"/R1.fastq -p "$threads" -k 32 2> "$qcOutDir"/mash.log

#identify species using genome database
echo "comparing to reference genome db"
mash screen -p "$threads" -w "$refDB" "$qcOutDir"/concatRawReads.fastq | sort -gr > "$qcOutDir"/mashscreen.genome.tsv

#identify if theres any plasmids present in the reads
echo "comparing to reference plasmid db"
mash screen -p "$threads" -w "$plasmidRefDB" "$qcOutDir"/concatRawReads.fastq | sort -gr > "$qcOutDir"/mashscreen.plasmid.tsv

source deactivate

#fastqc of the forward and reverse reads
echo "fastqc of reads"

source activate fastqc-0.11.7

fastqc "$R1" -o "$qcOutDir" --extract
fastqc "$R2" -o "$qcOutDir" --extract
rm -rf "$qcOutDir"/*.html
rm -rf "$qcOutDir"/*.zip

source deactivate

#kraken2 species classification
echo "kraken2 classfiication of reads"

source activate kraken2-2.0.7_beta

kraken2 --db "$k2db" --paired --threads "$threads" --report "$qcOutDir"/kraken2.genome.report --output "$qcOutDir"/kraken2.genome.classification "$R1" "$R2"

source deactivate

#calculate the total number of basepairs in the read
echo "calculating read stats"

source activate seqtk-1.3

R1bp=`seqtk fqchk $R1 | head -3 | tail -1 | cut -d$'\t' -f 2`
R2bp=`seqtk fqchk $R2 | head -3 | tail -1 | cut -d$'\t' -f 2`
totalbp=$((R1bp+R2bp))
echo $totalbp > "$qcOutDir"/totalbp

source deactivate

echo "cleaning up..."
rm -rf "$qcOutDir"/concatRawReads.fastq
rm -rf "$qcOutDir"/concatRawReads.fastq.msh
rm -rf "$qcOutDir"/R1.fastq
rm -rf "$qcOutDir"/R1.fastq.msh

echo "done step 1"

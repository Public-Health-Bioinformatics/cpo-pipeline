#!/bin/bash -e

#$ -V             # Pass environment variables to the job
#$ -N cpo_pipeline
#$ -cwd           # Use the current working dir
#$ -pe smp 8      # Parallel Environment (how many cores)
#$ -l h_vmem=11G  # Memory (RAM) allocation *per core*
#$ -e ./logs/$JOB_ID.err
#$ -o ./logs/$JOB_ID.log

# input parameters: 1=ID, 2=assemblyPath, 3=outputdir, 4=card.json 5=abricate_datadir 6=abricate_cpo_plasmid_db

ID="$1"
assemblyPath="$2"
outputDir="$3"
cardPath="$4"
abricate_datadir="$5"
abricate_carbapenemase_db="$6"
threads=8 #"$5"

echo "parameters: "
echo "ID: $ID"
echo "outputDir: $outputDir"
echo "card.json Path: $cardPath"
echo "threads: $threads"

typingDir="${outputDir}"/typing
resistanceDir="${outputDir}"/resistance

mkdir -p "${typingDir}"/"${ID}"/{"${ID}".mlst,"${ID}".recon}
mkdir -p "${resistanceDir}"/"${ID}"
mkdir -p "${outputDir}"/summary

echo "step3: mlst"

#step3: run mlst

source activate mlst-2.15.1

mlst "${assemblyPath}" > "${typingDir}/${ID}/${ID}.mlst/${ID}.mlst"

source deactivate

#step4: plasmid+amr
echo "step4: plasmid+amr prediction"

#find carbapenemases using custom cpo database.

source activate abricate-0.8.7

abricate --datadir "${abricate_datadir}" --db "${abricate_carbapenemase_db}" "${assemblyPath}" > "${resistanceDir}/${ID}/${ID}.cp"

source deactivate

#run rgi
cd "${resistanceDir}"/"${ID}"

source activate rgi-4.0.3

rgi load -i "${cardPath}" --local #--debug
rgi main -i "${assemblyPath}" -o "${ID}.rgi" -t contig -a BLAST -n "${threads}" --local --clean
rm -rf localDB

source deactivate

cd -

source activate mob_suite-1.4.5

#run mob typer
mob_recon --infile "${assemblyPath}" --outdir "${typingDir}/${ID}/${ID}.recon" --run_typer

source deactivate

echo "done step 4"

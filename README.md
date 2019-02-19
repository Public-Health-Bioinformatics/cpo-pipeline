[![Build Status](https://travis-ci.org/Public-Health-Bioinformatics/cpo-pipeline.svg?branch=master)](https://travis-ci.org/Public-Health-Bioinformatics/cpo-pipeline)

# cpo-pipeline

An analysis pipeline for the purpose of investigating [Carbapenemase-Producing Organisms](https://en.wikipedia.org/wiki/Carbapenem-resistant_enterobacteriaceae).

This codebase is derived from the [CPO_Prediction](https://github.com/imasianxd/CPO_Prediction) workflow, written by [Justin Jia](https://github.com/imasianxd) in collaboration with [William Hsiao](https://github.com/wwhsiao) and [Matt Croxen](https://github.com/mcroxen).

## Requirements

It is expected that several bioinformatics tools are provided as tool-specific conda environments. Each environment is activated in turn by calling:

```
source activate <tool-version>
```

...and deactivated by calling:

```
source deactivate
```

In order for this to work, a conda installation must be available on the user's `PATH`.

If these environments are not already available on your system, create them as follows:

```
conda create -n <tool-version> tool=version
```

For example:

```
conda create -n mash-2.0 mash=2.0
```

The following conda environments are required:

```
mash-2.0
seqtk-1.3
fastqc-0.11.8
shovill-1.0.1
quast-5.0.2
mlst-2.15.1
mob_suite-1.4.5
rgi-4.0.3
abricate-0.8.7
bwa-0.7.17
samtools-1.9
bcftools-1.9
freebayes-1.2
```

## Running the Pipeline

The pipeline can be run on a single sample as follows:

```
cpo-pipeline full --ID <SAMPLE_ID> --R1 <READ_1_FASTQ> --R2 <READ_2_FASTQ> --outdir <OUTPUT_DIRECTORY>
```

Pipeline output will be written to `OUTPUT_DIRECTORY/SAMPLE_ID` with the following sub-directories:

```
pre-assembly_qc
assembly
post-assembly_qc
reference
resistance
typing
plasmids
```



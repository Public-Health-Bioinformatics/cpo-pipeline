# CPO Pipeline Phase 1: Quality Control & Assembly

## Requirements

The `assembly/pipeline.py` script requires that the `pandas` library is installed for whichever `python` is invoked by `#!/usr/bin/env python`.

The `assembly/pipeline_qc.sh` script requires the following conda environments:

 - `mash-2.0`
 - `fastqc-0.11.7`
 - `seqtk-1.3`

The `assembly/pipeline_assembly.sh` script requires the following conda environments:

 - `shovill-1.0.1`
 - `quast-4.6.3`
 - `busco-3.0.2`

## Configuration

A configuration file located at `assembly/config.ini` Use this file to configure the install location of these scripts and the locations of several reference databases used. The config file is a simple key-value [ini](https://en.wikipedia.org/wiki/INI_file) format. eg:

```
[scripts]
script-path = /home/dfornika/code/cpo-pipeline/assembly

[databases]
mash-genomedb = /data/ref_databases/mash/refseq.genomes.k21s1000.msh
mash-plasmiddb = /data/ref_databases/mash/refseq.plasmid.k21s1000.msh
busco-db = /data/ref_databases/busco/enterobacteriales_odb9
```

Replace these paths with your system-specific paths.

## Running the Pipeline

### Setup

Prepare an analysis directory with three sub-directories: `input`, `output`, and `logs`.

The `input` directory should contain input `fastq[.gz]` files, along with a `sample_list.txt` file containing a simple list of sample IDs to include in the analysis.

```
Sample-001
Sample-002
Sample-003
```

```bash
.
├── input
│   ├── Sample-001_R1.fastq.gz
│   ├── Sample-001_R2.fastq.gz
│   ├── Sample-002_R1.fastq.gz
│   ├── Sample-002_R2.fastq.gz
│   ├── Sample-003_R1.fastq.gz
│   ├── Sample-003_R2.fastq.gz
│   ├── ...
│   └── sample_list.txt
├── logs
└── output
```

### Submit Analysis

```
assembly/pipeline_multisample.sh assembly/pipeline.qsub assembly/pipeline.py samples.txt input output "Expected species name"
```

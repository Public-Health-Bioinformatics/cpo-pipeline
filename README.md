[![Build Status](https://travis-ci.org/Public-Health-Bioinformatics/cpo-pipeline.svg?branch=master)](https://travis-ci.org/Public-Health-Bioinformatics/cpo-pipeline)

# cpo-pipeline

An analysis pipeline for the purpose of investigating [Carbapenemase-Producing Organisms](https://en.wikipedia.org/wiki/Carbapenem-resistant_enterobacteriaceae).

This codebase is derived from the [CPO_Prediction](https://github.com/imasianxd/CPO_Prediction) workflow, written by [Justin Jia](https://github.com/imasianxd) in collaboration with [William Hsiao](https://github.com/wwhsiao) and [Matt Croxen](https://github.com/mcroxen).

## Requirements

This pipeline uses [python bindings to DRMAA](https://github.com/pygridtools/drmaa-python), so it must be run on a DRMAA-compatible HPC compute cluster and the `DRMAA_LIBRARY_PATH` environment variable must be set. It is expected that several bioinformatics tools are provided as tool-specific conda environments which are available from the worker nodes of the cluster. Each environment is activated in turn by calling:

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

### Full Analysis, One Sample

The pipeline can be run on a single sample as follows:

```
cpo-pipeline full --ID <SAMPLE_ID> --R1 <READ_1_FASTQ> --R2 <READ_2_FASTQ> --outdir <OUTPUT_DIRECTORY>
```

### Full Analysis, Multiple Samples

The pipeline can be run on multiple samples in parallel as follows:

First, prepare a tab-delimited multi-sample input file with the following fields:

```
sample_id	path_to_read1_fastq_file	path_to_read2_fastq_file
```

You may optionally include a fourth field for the numeric NCBI taxonomy id for the expected species. If an expected species is provided, it will be used to estimate depth of sequence coverage based on the size of NCBI's representative genome for that species. If no expected species is specified, estimated coverage will be calculated based on the size of the reference genome that has the smallest mash distance from the sample.

To run the pipeline in multi-sample mode:

```
cpo-pipeline multi --input <MULTI_SAMPLE_INPUT_FILE> --outdir <OUTPUT_DIRECTORY>
```

If you have many samples to run and would like to limit the number of samples that are being analyzed concurrently, use the `-p` (or `--parallel`) flag:

```
cpo-pipeline multi --parallel 8 --input <MULTI_SAMPLE_INPUT_FILE> --outdir <OUTPUT_DIRECTORY>
```

### Partial Analyses

Each module of the pipeline can be run as an individual unit. Each unit accepts a `-h` or `--help` flag that will provide usage information. To see all available sub-commands:

```
$ cpo-pipeline --help
usage: cpo-pipeline [-h] [-v] {full,multi,assembly,plasmids,typing,resistance,tree} ...

positional arguments:
  {full,multi,assembly,plasmids,typing,resistance,tree}

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

To see usage for an individual sub-command:

```
$ cpo-pipeline plasmids --help
usage: cpo-pipeline plasmids [-h] -i SAMPLE_ID -1 READS1_FASTQ -2 READS2_FASTQ [-o OUTDIR]
                             [--mash-refseq-plasmid-db MASH_REFSEQ_PLASMID_DB]
                             [--mash-custom-plasmid-db MASH_CUSTOM_PLASMID_DB] [-c CONFIG_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -i SAMPLE_ID, --ID SAMPLE_ID
                        identifier of the isolate
  -1 READS1_FASTQ, --R1 READS1_FASTQ
                        absolute file path forward read (R1)
  -2 READS2_FASTQ, --R2 READS2_FASTQ
                        absolute file path to reverse read (R2)
  -o OUTDIR, --outdir OUTDIR
                        absolute path to output folder
  --mash-refseq-plasmid-db MASH_REFSEQ_PLASMID_DB
                        absolute path to refseq mash reference database
  --mash-custom-plasmid-db MASH_CUSTOM_PLASMID_DB
                        absolute file path to directory of custom plasmid mash sketches
  -c CONFIG_FILE, --config CONFIG_FILE
                        Config File
```

## Pipeline Output

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

Final outputs are collected for each sample in the files: `final_output.tsv` and `final_plasmid.tsv`

## Logging

The pipeline emits json-formatted structured logs using the `structlog` library. Each log entry represents an event that occurrs while the pipeline is running.

| event type                     |
|--------------------------------|
| `analysis_started`             |
| `job_submitted`                |
| `job_completed`                |
| `parsed_result_file`           |
| `download_failed`              |
| `result_parsing_failed`        |
| `failed_quality_control_check` |
| `collected_final_outputs`      |


All events include the following fields:

| field                    | type      | description                                      |
|--------------------------|-----------|--------------------------------------------------|
| `pipeline_version`       | `string`  | semantic versioning (eg. `0.1.0`)                |
| `level`                  | `string`  | log level (`info` or `error`)                    |
| `analysis_id`            | `uuid4`   | unique id for each sample, for each analysis run |
| `timestamp`              | `string`  | ISO-8601 format, microsecond precision, UTC zone YYYY-MM-DDTHH:mm:ss.xxxxxx+00:00 |
| `event`                  | `string`  | one of the events listed in the table above      |

A typical log entry looks like this (although it won't be pretty-printed or formatted during an actual pipeline run):

```json
{
  "analysis_id": "883290b2-0f42-476c-94d7-e18cdddb8784",
  "sample_id": "BC18-Cfr003",
  "pipeline_version": "0.1.0",
  "timestamp": "2019-04-12T23:55:00.790729+00:00",
  "event": "analysis_started",
  "level": "info"
}
```

The logs may be a bit difficult for a person to visually parse while the pipeline is running, but they're designed to be machine-readable. They can also potentially be fed into a system like the [ELK stack](https://www.elastic.co/elk-stack) or queried by tools like [Apache Drill](https://drill.apache.org/docs/querying-json-files/).

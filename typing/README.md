# CPO Pipeline Phase 2: MLST Typing, Resistance Gene Identification, Plasmid Analysis

## Requirements

The `typing/pipeline.py` script requires that the `pandas` library is installed for whichever `python` is invoked by `#!/usr/bin/env python`.

The `typing/pipeline_typing.qsub` script requires the following conda environments:

 - `mlst-2.15.1`
 - `abricate-0.8.7`
 - `rgi-4.0.3`
 - `mob_suite-1.4.5`

## Configuration

A configuration file located at `typing/config.ini` Use this file to configure the install location of these scripts and the locations of several reference databases used. The config file is a simple key-value [ini](https://en.wikipedia.org/wiki/INI_file) format. eg:

```
[scripts]
script-path = /home/dfornika/code/cpo-pipeline/typing

[databases]
card = /data/ref_databases/card/card-2.0.3.json
abricate-datadir = /projects/carbapenemase_producing_organisms_surveillance/databases/abricate
abricate-cpo-plasmid-db = carbapenemases
mlst-scheme-map = /home/dfornika/code/cpo-pipeline/typing/scheme_species_map.tab
```

Note: An abricate database containing carbapenemase gene sequences must be prepared before running the pipeline. (todo: provide details on how to build that db)

### Setup

This module uses outputs from [Module 1 (Assembly & QC)](https://github.com/Public-Health-Bioinformatics/cpo-pipeline/tree/master/assembly). That module must be run first before running this one.

As with Module 1, this module requires a `sample_list.txt` file containing a simple list of sample IDs to include in the analysis:

```
Sample-001
Sample-002
Sample-003
```

### Submit Analysis

1. Run `mlst` typing, `abricate` carbapenemase gene identification, `rgi` resistance gene identification and `mob_suite` mobile element (plasmid analysis):

```
typing/pipeline_multisample.sh typing/pipeline_typing.qsub sample_list.txt module_1_output_dir module_2_output_dir
```

2. 

```
typing/pipeline.py
```
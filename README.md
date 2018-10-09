# cpo-pipeline

An analysis pipeline for the purpose of investigating [Carbapenemase-Producing Organisms](https://en.wikipedia.org/wiki/Carbapenem-resistant_enterobacteriaceae).

# Module 1: Quality Control & Assembly

## Requirements

The `assembly/pipeline.py` script requires that the `pandas` library is installed for whichever `python` is invoked by `#!/usr/bin/env python`.

The `assembly/pipeline_qc.sh` script requires the following conda environments:

 - `mash-2.0`
 - `fastqc-0.11.7`
 - `kraken2-2.0.7_beta`
 - `seqtk-1.3`

The `assembly/pipeline_assembly.sh` script requires the following conda environments:

 - `shovill-1.0.1`
 - `quast-4.6.3`
 - `busco-3.0.2`

The `assembly/pipeline_assembly_contaminant.sh` script requires the following conda environments:

 - `bbmap-38.22`
 - `shovill-1.0.1`
 - `quast-4.6.3`
 - `busco-3.0.2`

# Module 2: MLST Typing, Resistance Gene Identification, Plasmid Analysis and Dendrogram

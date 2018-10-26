# CPO Pipeline Phase 3: Dendrogram (Tree)

## Requirements

The `tree/pipeline.py` script requires that the `pandas` library is installed for whichever `python` is invoked by `#!/usr/bin/env python`.

The `tree/pipeline_tree.sh` script requires the following conda environments:

 - `snippy-4.2.3`

## Configuration

A configuration file located at `tree/config.ini` Use this file to configure the install location of these scripts and the locations of several reference databases used. The config file is a simple key-value [ini](https://en.wikipedia.org/wiki/INI_file) format. eg:

```
[scripts]
script-path = /home/dfornika/code/cpo-pipeline/tree
```

### Setup

This module uses outputs from [Phase 1 (Assembly & QC)](https://github.com/Public-Health-Bioinformatics/cpo-pipeline/tree/master/assembly) and [Phase 2 (Typing & Resistance)](https://github.com/Public-Health-Bioinformatics/cpo-pipeline/tree/master/typing). Those pipelines must be run first before running this one.

As with Phases 1 and 2, this module requires a `sample_list.txt` file containing a simple list of sample IDs to include in the analysis:

```
Sample-001
Sample-002
Sample-003
```

### Submit Analysis

TBD
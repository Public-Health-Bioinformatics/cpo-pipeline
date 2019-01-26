# -*- coding: utf-8 -*-
"""cpo-pipeline.assembly.parsers.input_parsers

This module provides functions for parsing input files provided to
the Assembly phase of the cpo-pipeline.
"""

import csv

def parse_estimated_genome_sizes(input_path):
    estimated_genome_sizes = []
    with open(input_path) as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        for row in reader:
            row['estimated_genome_size'] = int(row['estimated_genome_size'])
            estimated_genome_sizes.append(row)

    return estimated_genome_sizes

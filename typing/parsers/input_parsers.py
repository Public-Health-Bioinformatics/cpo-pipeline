# -*- coding: utf-8 -*-
"""cpo-pipeline.typing.parsers.input_parsers

This module provides functions for parsing input files provided to
the Typing phase of the cpo-pipeline.
"""

import csv

def parse_scheme_species_map(path_to_scheme_species_map):
    """
    Args:
        path_to_scheme_species_map (str): Path to the scheme-species map file.

    Returns:
        list of dict: Parsed mlst report.
        For example:
        [
            {
                'scheme_id': 'abaumannii',
                'genus': 'Acinetobacter',
                'species': 'Acinetobacter baumannii',
            },
            {
                'scheme_id': 'abaumannii_2',
                'genus': 'Acinetobacter',
                'species': 'Acinetobacter baumannii',
            },
            }
                'scheme_id': 'achromobacter',
                'genus': 'Achromobacter'
            },
            ...
        ]
        Note: some schemes only correspond to a genus, not a species.
    """
    scheme_species_map = []
    with open(path_to_scheme_species_map) as scheme_species_map_file:
        reader = csv.reader(scheme_species_map_file, delimiter='\t')
        next(reader) # skip header
        for row in reader:
            scheme_species_map_record = {}
            scheme_species_map_record['scheme_id'] = row[0]
            scheme_species_map_record['genus'] = row[1]
            if row[2] != '':
                scheme_species_map_record['species'] = row[1] + " " + row[2]
            scheme_species_map.append(scheme_species_map_record)
    # TODO: Validate that each entry has a unique scheme_id
    return scheme_species_map

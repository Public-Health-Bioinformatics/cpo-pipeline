import os
import sys
import unittest
import json
import pprint

# add the '../../typing' directory to the 'sys.path'
TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(TEST_DIR_PATH))

from parsers import result_parsers
from parsers import input_parsers

class SchemeSpeciesMapParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_input/scheme_species_map.tsv')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_input/parsed_scheme_species_map.json')) as scheme_species_map_json_file:
            self.parsed_scheme_species_map_json = json.load(scheme_species_map_json_file)
            scheme_species_map_json_file.close()

    def test_parse_mlst_result(self):
        parsed_result = input_parsers.parse_scheme_species_map(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_scheme_species_map_json)
        for parsed_scheme_species_map_record, parsed_scheme_species_map_json_record in paired_results:
            self.assertDictEqual(parsed_scheme_species_map_record, parsed_scheme_species_map_json_record)

            
class MlstResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/typing/SAMPLE-ID/SAMPLE-ID.mlst/SAMPLE-ID.mlst')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_mlst_report.json')) as mlst_report_json:
            self.parsed_mlst_report_json = json.load(mlst_report_json)
            mlst_report_json.close()

    def test_parse_mlst_result(self):
        parsed_result = result_parsers.parse_mlst_result(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_mlst_report_json)
        for parsed_mlst_report_record, parsed_mlst_json_record in paired_results:
            self.assertDictEqual(parsed_mlst_report_record, parsed_mlst_json_record)
        
        
if __name__ == '__main__':
    unittest.main()

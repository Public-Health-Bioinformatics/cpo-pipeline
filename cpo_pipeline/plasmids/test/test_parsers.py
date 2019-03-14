import os
import sys
import unittest
import json
import pprint

TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(TEST_DIR_PATH))

import parsers

class PlasmidDataParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_input/NZ_CP003997.dat')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_input/NZ_CP003997.dat.json')) as plasmid_data_json:
            self.parsed_plasmid_data_json = json.load(plasmid_data_json)
            plasmid_data_json.close()

    def test_custom_plasmid_db_dat_parser(self):
        parsed_result = parsers.custom_plasmid_db_dat_parser(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_plasmid_data_json)
        for parsed_plasmid_data_record, parsed_plasmid_data_json_record in paired_results:
            self.assertDictEqual(parsed_plasmid_data_record, parsed_plasmid_data_json_record)


if __name__ == '__main__':
    unittest.main()

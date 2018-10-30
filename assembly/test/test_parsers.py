import os
import sys
import unittest
import json

# add the '../../assembly' directory to the 'sys.path'

TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(TEST_DIR_PATH))

from parsers import result_parsers

class KrakenResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/kraken_report.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_kraken_report.json')) as kraken_report_json:
            self.parsed_kraken_json = json.load(kraken_report_json)
            kraken_report_json.close()
        
    def test_parse_kraken_result(self):
        parsed_result = result_parsers.parse_kraken_result(self.test_data_path)
        self.assertDictEqual(parsed_result, self.parsed_kraken_json)
    
if __name__ == '__main__':
    unittest.main()

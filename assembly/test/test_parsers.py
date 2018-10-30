import os
import sys
import unittest
import json

# add the '../../assembly' directory to the 'sys.path'
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers import result_parsers

class KrakenResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/kraken_report.txt')
        with open('data/parsed_kraken_report.json') as kraken_report_json:
            self.parsed_kraken_json = json.load(kraken_report_json)
            kraken_report_json.close()
        
    def test_parse_kraken_result(self):
        parsed_result = result_parsers.parse_kraken_result(self.test_data_path)
        self.assertDictEqual(parsed_result, self.parsed_kraken_json)
    
if __name__ == '__main__':
    unittest.main()

import os
import sys
import unittest
import json
from collections import OrderedDict

# add the '../../resistance' directory to the 'sys.path'
TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(TEST_DIR_PATH))

from parsers import result_parsers

class AbricateResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/SAMPLE-ID/resistance/abricate/abricate.tsv')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_abricate_report.json')) as abricate_report_json_file:
            self.parsed_abricate_report_json = json.load(abricate_report_json_file)
            abricate_report_json_file.close()

    def test_parse_abricate_result(self):
        parsed_result = result_parsers.parse_abricate_result(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_abricate_report_json)
        for parsed_abricate_report_record, parsed_abricate_json_record in paired_results:
            self.assertDictEqual(parsed_abricate_report_record, parsed_abricate_json_record)

class RgiTxtReportParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/SAMPLE-ID/resistance/rgi/rgi.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_rgi_txt_report.json')) as rgi_txt_report_json_file:
            self.parsed_rgi_txt_report_json = json.load(rgi_txt_report_json_file, object_pairs_hook=OrderedDict)
            rgi_txt_report_json_file.close()

    def test_parse_rgi_result_txt(self):
        parsed_result = result_parsers.parse_rgi_result_txt(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_rgi_txt_report_json)
        for parsed_rgi_txt_report_record, parsed_rgi_txt_report_json_record in paired_results:
            self.assertDictEqual(parsed_rgi_txt_report_record, parsed_rgi_txt_report_json_record)

if __name__ == '__main__':
    unittest.main()

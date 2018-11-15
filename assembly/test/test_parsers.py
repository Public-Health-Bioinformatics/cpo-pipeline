import os
import sys
import unittest
import json
import pprint

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
        paired_results = zip(parsed_result, self.parsed_kraken_json)
        for parsed_kraken_report_record, parsed_kraken_json_record in paired_results:
            self.assertDictEqual(parsed_kraken_report_record, parsed_kraken_json_record)


class FastqcResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/fastqc_summary.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_fastqc_summary.json')) as fastqc_summary_json:
            self.parsed_fastqc_summary_json = json.load(fastqc_summary_json)
            fastqc_summary_json.close()

    def test_parse_fastqc_result(self):
        parsed_result = result_parsers.parse_fastqc_result(self.test_data_path)
        self.assertDictEqual(parsed_result, self.parsed_fastqc_summary_json)

class MashResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/mashscreen.genome.tsv')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_mash_screen_genome_report.json')) as mash_screen_genome_json:
            self.parsed_mash_screen_genome_json = json.load(mash_screen_genome_json)
            mash_screen_genome_json.close()

    def test_parse_mash_result(self):    
        parsed_result = result_parsers.parse_mash_result(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_mash_screen_genome_json)
        for parsed_mash_report_record, parsed_mash_json_record in paired_results:
            self.assertDictEqual(parsed_mash_report_record, parsed_mash_json_record)

class MashLogParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/mash.log')
        self.test_totalbp_path = os.path.join(TEST_DIR_PATH, 'data/totalbp')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_mash_log.json')) as mash_log_json:
            self.parsed_mash_log_json = json.load(mash_log_json)
            mash_log_json.close()

    def test_parse_read_stats(self):
        parsed_result = result_parsers.parse_read_stats(self.test_data_path, self.test_totalbp_path)
        self.assertDictEqual(parsed_result, self.parsed_mash_log_json)

if __name__ == '__main__':
    unittest.main()

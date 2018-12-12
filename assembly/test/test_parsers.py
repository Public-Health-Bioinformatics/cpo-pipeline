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
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/kraken2.genome.report')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_kraken_report.json')) as kraken_report_json:
            self.parsed_kraken_json = json.load(kraken_report_json)
            kraken_report_json.close()
        
    def test_parse_kraken_result(self):
        parsed_result = result_parsers.parse_kraken_result(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_kraken_json)
        for parsed_kraken_report_record, parsed_kraken_json_record in paired_results:
            self.assertDictEqual(parsed_kraken_report_record, parsed_kraken_json_record)


class FastqcResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/SAMPLE-ID_R1_fastqc/summary.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_fastqc_summary.json')) as fastqc_summary_json:
            self.parsed_fastqc_summary_json = json.load(fastqc_summary_json)
            fastqc_summary_json.close()

    def test_parse_fastqc_result(self):
        parsed_result = result_parsers.parse_fastqc_result(self.test_data_path)
        self.assertDictEqual(parsed_result, self.parsed_fastqc_summary_json)

class MashResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/mashscreen.genome.tsv')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_mash_screen_genome_report.json')) as mash_screen_genome_json:
            self.parsed_mash_screen_genome_json = json.load(mash_screen_genome_json)
            mash_screen_genome_json.close()

    def test_parse_mash_result(self):    
        parsed_result = result_parsers.parse_mash_result(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_mash_screen_genome_json)
        for parsed_mash_report_record, parsed_mash_json_record in paired_results:
            self.assertDictEqual(parsed_mash_report_record, parsed_mash_json_record)

class MashLogParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/mash.log')
        self.test_totalbp_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/totalbp')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_mash_log.json')) as mash_log_json:
            self.parsed_mash_log_json = json.load(mash_log_json)
            mash_log_json.close()

    def test_parse_read_stats(self):
        parsed_result = result_parsers.parse_read_stats(self.test_data_path, self.test_totalbp_path)
        self.assertDictEqual(parsed_result, self.parsed_mash_log_json)

class TotalBpParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/totalbp')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_totalbp.json')) as totalbp_json:
            self.parsed_totalbp_json = json.load(totalbp_json)
            totalbp_json.close()
    
    def test_parse_total_bp(self):
        parsed_result = result_parsers.parse_total_bp(self.test_data_path)
        self.assertEqual(parsed_result, self.parsed_totalbp_json['total_bp'])

class ReferenceGenomeStatsParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/qcResult/SAMPLE-ID/GCF_001022155.1_ASM102215v1_assembly_stats.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_reference_genome_stats.json')) as reference_genome_stats_json:
            self.parsed_reference_genome_stats_json = json.load(reference_genome_stats_json)
            reference_genome_stats_json.close()
    
    def test_parse_reference_genome_stats(self):
        parsed_result = result_parsers.parse_reference_genome_stats(self.test_data_path)
        self.assertEqual(parsed_result, self.parsed_reference_genome_stats_json['reference_genome_total_length'])

class BuscoResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/assembly_qc/SAMPLE-ID/SAMPLE-ID.busco/short_summary_SAMPLE-ID.busco.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_busco_summary.json')) as busco_summary_json:
            self.parsed_busco_summary_json = json.load(busco_summary_json)
            busco_summary_json.close()

    def test_parse_busco_result(self):
        parsed_result = result_parsers.parse_busco_result(self.test_data_path)
        self.assertDictEqual(parsed_result, self.parsed_busco_summary_json)

class QuastResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/assembly_qc/SAMPLE-ID/SAMPLE-ID.quast/report.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_quast_report.json')) as quast_report_json:
            self.parsed_quast_report_json = json.load(quast_report_json)
            quast_report_json.close()
    
    def test_parse_quast_result(self):
        parsed_result = result_parsers.parse_quast_result(self.test_data_path)
        self.assertDictEqual(parsed_result, self.parsed_quast_report_json)


if __name__ == '__main__':
    unittest.main()

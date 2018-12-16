import os
import sys
import unittest
import json
import pprint
from collections import OrderedDict

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
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/SAMPLE-ID/typing/mlst/mlst.tsv')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_mlst_report.json')) as mlst_report_json:
            self.parsed_mlst_report_json = json.load(mlst_report_json)
            mlst_report_json.close()

    def test_parse_mlst_result(self):
        parsed_result = result_parsers.parse_mlst_result(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_mlst_report_json)
        for parsed_mlst_report_record, parsed_mlst_json_record in paired_results:
            self.assertDictEqual(parsed_mlst_report_record, parsed_mlst_json_record)

class MobReconContigReportParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/SAMPLE-ID/typing/mob_recon/contig_report.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_mob_recon_contig_report.json')) as mob_recon_contig_report_json_file:
            self.parsed_mob_recon_contig_report_json = json.load(mob_recon_contig_report_json_file)
            mob_recon_contig_report_json_file.close()

    def test_parse_mob_recon_contig_report(self):
        parsed_result = result_parsers.parse_mob_recon_contig_report(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_mob_recon_contig_report_json)
        for parsed_mob_recon_contig_report_record, parsed_mob_recon_contig_report_json_record in paired_results:
            self.assertDictEqual(parsed_mob_recon_contig_report_record, parsed_mob_recon_contig_report_json_record)

class MobReconMobTyperAggregateReportParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/SAMPLE-ID/typing/mob_recon/mobtyper_aggregate_report.txt')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_mob_recon_mobtyper_aggregate_report.json')) as mob_recon_mobtyper_aggregate_report_json_file:
            self.parsed_mob_recon_mobtyper_aggregate_report_json = json.load(mob_recon_mobtyper_aggregate_report_json_file, object_pairs_hook=OrderedDict)
            mob_recon_mobtyper_aggregate_report_json_file.close()

    def test_parse_mob_recon_mobtyper_aggregate_report(self):
        parsed_result = result_parsers.parse_mob_recon_mobtyper_aggregate_report(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_mob_recon_mobtyper_aggregate_report_json)
        for parsed_mob_recon_mobtyper_aggregate_report_record, parsed_mob_recon_mobtyper_aggregate_report_json_record in paired_results:
            self.assertDictEqual(parsed_mob_recon_mobtyper_aggregate_report_record, parsed_mob_recon_mobtyper_aggregate_report_json_record)

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

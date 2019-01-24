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

class WorkflowResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(TEST_DIR_PATH, 'data/pipeline_output/workflow_results.tsv')
        with open(os.path.join(TEST_DIR_PATH, 'data/parsed_results/parsed_workflow_results.json')) as workflow_results_json_file:
            self.parsed_workflow_results_json = json.load(workflow_results_json_file)
            workflow_results_json_file.close()

    def test_parse_workflow_results(self):
        parsed_result = result_parsers.parse_workflow_results(self.test_data_path)
        paired_results = zip(parsed_result, self.parsed_workflow_results_json)
        for parsed_workflow_results_record, parsed_workflow_results_json_record in paired_results:
            self.assertDictEqual(parsed_workflow_results_record, parsed_workflow_results_json_record)

if __name__ == '__main__':
    unittest.main()

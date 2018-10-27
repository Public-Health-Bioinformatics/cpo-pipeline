import os
import sys
import unittest

print( os.path.dirname(os.path.dirname(os.path.abspath(__file__))) )

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from parsers import result_parsers

class KrakenResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(os.path.dirname(__file__), '/data/kraken_report.txt')

    def test_parse_kraken_result(self):
        parsed_result = result_parsers.parse_kraken_result(self.test_data_path)
        self.assertTrue(False) # TBD
    

import os
import sys
import unittest

# add the '../../assembly' directory to the 'sys.path'
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers import result_parsers

class KrakenResultParserTest(unittest.TestCase):
    def setUp(self):
        self.test_data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/kraken_report.txt')
        print(self.test_data_path)
        
    def test_parse_kraken_result(self):
        parsed_result = result_parsers.parse_kraken_result(self.test_data_path)
        self.assertTrue(False) # TBD
    
if __name__ == '__main__':
    unittest.main()

import os
import sys
import unittest
import json
import pprint

# add the '../../typing' directory to the 'sys.path'
TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(TEST_DIR_PATH))

from parsers import result_parsers

if __name__ == '__main__':
    unittest.main()

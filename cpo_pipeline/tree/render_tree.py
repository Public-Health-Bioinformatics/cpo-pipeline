#!/usr/bin/env python

import argparse
import os
import sys
os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree, TreeStyle, NodeStyle

def main(args):
    tree = Tree(args.input_file, format=1)
    tree.set_outgroup(tree&"Reference")
    tree_style = TreeStyle()
    tree_style.show_leaf_name = True
    tree.render(args.output_file, tree_style=tree_style)

if __name__ == '__main__':
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name, description='')
    parser.add_argument("-i", "--input", dest="input_file",
                        help="tree file", required=True)
    parser.add_argument("-o", "--output", dest="output_file",
                        help="tree file", required=True)
    args = parser.parse_args()
    main(args)

#!/usr/bin/env python3

from Bio.Phylo import read as tree_read
import sys

def main():
    tree = tree_read(sys.argv[1], 'newick')
    tree.root_at_midpoint()
    for leaf in tree.get_terminals():
        print(leaf.name, file=sys.stdout)

if __name__ == "__main__":
    main()

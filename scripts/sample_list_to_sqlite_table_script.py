#!/usr/bin/env python3

import sys

if __name__ == "__main__":
    sample_list_path = sys.argv[1]
    k = int(sys.argv[2])
    with open(sample_list_path) as f:
        sample_columns = ', '.join([line.split()[0] + ' INT' for line in f])
    output = f"CREATE TABLE count_ (kmer VARCHAR({k}) PRIMARY KEY, {sample_columns});"
    print(output)

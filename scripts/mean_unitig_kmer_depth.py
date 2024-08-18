#!/usr/bin/env python3

import sys
import numpy as np
import strainzip as sz
from lib.util import info
from tqdm import tqdm
from collections import defaultdict

if __name__ == "__main__":
    k_size = int(sys.argv[1])
    n_samples = int(sys.argv[2])
    fasta_inpath = sys.argv[3]

    info("Loading unitigs")
    kmer_to_unitig_id = {}
    unitig_kmer_length = {}
    with open(fasta_inpath) as f:
        for header, sequence in tqdm(sz.io.iter_linked_fasta_entries(f)):
            unitig_id_string, *_ = sz.io.generic_header_tokenizer(header)
            unitig_id = unitig_id_string[1:]
            unitig_kmer_length[unitig_id] = len(sequence) - k_size + 1
            for kmer in sz.sequence.iter_kmers(sequence, k=k_size):
                kmer_to_unitig_id[kmer] = unitig_id

    info("Accumulating kmer counts for each unitig.")
    unitig_kmer_count = defaultdict(lambda: np.zeros(n_samples, dtype="int"))
    for line in tqdm(sys.stdin):
        kmer, *counts = line.split("\t")
        counts = [int(token.strip()) for token in counts]
        if kmer not in kmer_to_unitig_id:  # Not canonical
            kmer = sz.io.reverse_complement(kmer)
        unitig_kmer_count[kmer_to_unitig_id[kmer]] += counts

    info("Calculating unitig mean depth and writing output.")
    for unitig_id in unitig_kmer_length:
        _unitig_depth = unitig_kmer_count[unitig_id] / unitig_kmer_length[unitig_id]
        print(unitig_id, *_unitig_depth, sep="\t")

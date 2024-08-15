#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import strainzip as sz
from lib.util import info
from tqdm import tqdm
from collections import defaultdict

if __name__ == "__main__":
    k = int(sys.argv[1])
    sample_names = sys.argv[2].split(",")
    fasta_inpath = sys.argv[3]
    outpath = sys.argv[4]

    n_samples = len(sample_names)

    info("Loading unitigs")
    kmer_to_unitig_id = {}
    unitig_kmer_length = {}
    with open(fasta_inpath) as f:
        for header, sequence in tqdm(sz.io.iter_linked_fasta_entries(f)):
            unitig_id_string, *_ = sz.io.generic_header_tokenizer(header)
            unitig_id = unitig_id_string[1:]
            unitig_kmer_length[unitig_id] = len(sequence) - k + 1
            for kmer in sz.sequence.iter_kmers(sequence, k=k):
                kmer_to_unitig_id[kmer] = unitig_id

    unitig_kmer_length = pd.Series(unitig_kmer_length)

    info("Accumulating kmer counts for each unitig.")
    unitig_kmer_count = defaultdict(lambda: np.zeros(n_samples, dtype='int'))
    for line in tqdm(sys.stdin):
        kmer, *counts = line.split("\t")
        counts = [int(token.strip()) for token in counts]
        if kmer not in kmer_to_unitig_id:  # Not canonical
            kmer = sz.io.reverse_complement(kmer)
        unitig_kmer_count[kmer_to_unitig_id[kmer]] += counts

    info("Calculating mean depth.")
    unitig_depth = (
        pd.DataFrame(
            unitig_kmer_count.values(),
            index=unitig_kmer_count.keys(),
            columns=sample_names,
        )
        .divide(unitig_kmer_length, axis=0)
        .fillna(0)
        .rename(index=int)  # Unitigs are ints
        .rename_axis(index="unitig", columns="sample")
        .stack()
        .to_xarray()
    )

    info("Writing output.")
    unitig_depth.to_netcdf(outpath)

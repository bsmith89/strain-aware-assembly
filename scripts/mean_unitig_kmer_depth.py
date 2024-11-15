#!/usr/bin/env python3

import sys
import numpy as np
import strainzip as sz
from lib.util import info
from tqdm import tqdm
from collections import defaultdict
import xarray as xr

if __name__ == "__main__":
    k_size = int(sys.argv[1])
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
            unitig_kmer_length[unitig_id] = len(sequence) - k_size + 1
            for kmer in sz.sequence.iter_kmers(sequence, k=k_size):
                kmer_to_unitig_id[kmer] = unitig_id

    info("Accumulating kmer counts for each unitig.")
    unitig_kmer_count = defaultdict(lambda: np.zeros(n_samples, dtype="int"))
    # TODO: Consider having just one working "counts" buffer, that gets
    # over-written for each kmer, instead of allocating every single time.
    for line in tqdm(sys.stdin):
        kmer, *counts = line.split("\t")
        # TODO: Convert this list comprehension to a np.array(counts).astype('int')
        counts = [int(token.strip()) for token in counts]
        if kmer not in kmer_to_unitig_id:  # Not canonical
            kmer = sz.io.reverse_complement(kmer)
        unitig_kmer_count[kmer_to_unitig_id[kmer]] += counts

    info("Calculating unitig mean depth.")
    unitig_kmer_depth = []
    for unitig_id in tqdm(unitig_kmer_length):
        # NOTE: Appending to the matrix and then delete from the counts
        # to keep memory usage in check.
        unitig_kmer_depth.append(
            unitig_kmer_count[unitig_id] / unitig_kmer_length[unitig_id]
        )
        del unitig_kmer_count[unitig_id]

    unitig_kmer_depth = xr.DataArray(
        np.stack(unitig_kmer_depth),
        coords=dict(unitig=np.array(list(unitig_kmer_length.keys())).astype(int), sample=sample_names),
    )

    info("Writing output.")
    unitig_kmer_depth.to_netcdf(outpath)

#!/usr/bin/env python3

import sys
from math import ceil
from collections import defaultdict


if __name__ == "__main__":
    tile_len = int(sys.argv[1])
    overlap_len = int(sys.argv[2])
    inpath = sys.argv[3]

    raw_seqs = defaultdict(lambda: [])  # Sequences stored as list of lines
    with open(inpath) as f:
        for line in f:
            if line.startswith(">"):
                seqid = line.split()[0][1:]  # Sequence ID is the first word on the header line.
            else:
                raw_seqs[seqid].append(line.strip())

    tile_offset = tile_len - overlap_len
    for seqid, seq_lines in raw_seqs.items():
        seq = "".join(seq_lines)
        seq_len = len(seq)
        n_tiles = ceil((seq_len - tile_len) / tile_offset)
        stop = 0
        for i in range(n_tiles):
            start = i * tile_offset
            stop = i * tile_offset + tile_len
            print(f">{seqid}[{start}-{stop}]+")
            print(seq[start:stop])
        if stop != seq_len:
            stop = seq_len
            start = stop - tile_len
            print(f">{seqid}[{start}-{stop}]")
            print(seq[start:stop])
        # print(seq_len, n_tiles, file=sys.stderr)

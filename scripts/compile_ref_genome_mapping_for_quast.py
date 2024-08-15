#!/usr/bin/env python3

import pandas as pd
import sys
from tqdm import tqdm

def collect_fasta_entry_headers(handle):
    for line in handle:
        if line.startswith('>'):
            contig_name = line[1:].split()[0]
            yield contig_name

if __name__ == "__main__":
    genome_inpaths = dict(arg.split('=') for arg in sys.argv[1:])
    contigs = []
    for genome in tqdm(genome_inpaths):
        with open(genome_inpaths[genome]) as f:
            d0 = pd.Series(collect_fasta_entry_headers(f)).to_frame('contig')
            d1 = d0.assign(genome=genome)[['contig', 'genome']]
        contigs.append(d1)
    contigs = pd.concat(contigs).assign(quast_contig=lambda x: x.genome + '_' + x.contig)
    contigs[['quast_contig', 'genome']].to_csv(sys.stdout, sep='\t', index=False)


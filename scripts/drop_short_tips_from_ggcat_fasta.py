#!/usr/bin/env python3

import sys
import strainzip as sz
from tqdm import tqdm


if __name__ == "__main__":
    k = int(sys.argv[1])
    length_mult = float(sys.argv[2])
    inpath = sys.argv[3]
    outpath = sys.argv[4]

    min_length = k * length_mult
    header_tokenizer = sz.io.ggcat_header_tokenizer

    with open(inpath) as in_handle, open(outpath, 'w') as out_handle:
        for header, sequence in tqdm(sz.io.iter_linked_fasta_entries(in_handle)):
            # Parse the header.
            unitig_id_string, length_string, edge_strings_list = header_tokenizer(header)
            unitig_id, length, edge_list = sz.io.parse_linked_fasta_entry_header(
                unitig_id_string, length_string, edge_strings_list, k, header_tokenizer
            )

            # Determine if it's a tip.
            is_tip = False
            if len(edge_list) < 2:
                is_tip = True
            else:
                edge_directions = set([e[0][-1] for e in edge_list])
                if not '+' in edge_directions:
                    is_tip = True
                if not '-' in edge_directions:
                    is_tip = True

            # Write out any unitig that's either not a tip or sufficiently long.
            if (length >= min_length) or (not is_tip):
                print(f'>{unitig_id}\n{sequence}', file=out_handle)

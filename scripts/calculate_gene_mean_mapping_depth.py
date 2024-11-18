#!/usr/bin/env python3
# USAGE: {input.script} {input.gff} {input.pos_depth} {output}
import sys
import pandas as pd
from itertools import groupby
from tqdm import tqdm


def parse_gff_to_genes(text_stream):
    for line in text_stream:
        if line.startswith("#"):
            continue
        (
            seqid,
            source,
            _type,
            left_p1,
            right,
            score,
            strand,
            phase,
            attributes,
        ) = line.strip().split("\t")
        yield seqid, int(left_p1) - 1, int(right), strand


def group_and_tabularize_genes_by_seqid(gene_record_stream):
    for seqid, chunk in groupby(gene_record_stream, key=lambda record: record[0]):
        yield seqid, pd.DataFrame(chunk, columns=["seqid", "left", "right", "strand"])


def parse_depths(text_stream):
    for line in text_stream:
        seqid, pos, depth = line.strip().split("\t")
        yield seqid, int(pos), int(depth)


def group_and_tabularize_depths_by_seqid(depth_record_stream):
    for seqid, chunk in groupby(depth_record_stream, key=lambda record: record[0]):
        yield seqid, pd.DataFrame(chunk, columns=["seqid", "pos", "depth"])


def join_chunk_iters(gene_group_stream, depth_group_stream):
    "Assumes that both streams are sorted in ascending order by seqid."
    gene_seqid, gene_group = next(gene_group_stream)
    depth_seqid, depth_group = next(depth_group_stream)
    while True:
        if gene_seqid < depth_seqid:
            try:
                gene_seqid, gene_group = next(gene_group_stream)
            except StopIteration:
                return
        elif gene_seqid > depth_seqid:
            try:
                depth_seqid, depth_group = next(depth_group_stream)
            except StopIteration:
                return
        else:
            assert gene_seqid == depth_seqid
            yield gene_seqid, gene_group, depth_group
            try:
                gene_seqid, gene_group = next(gene_group_stream)
                depth_seqid, depth_group = next(depth_group_stream)
            except StopIteration:
                return


def calculate_mean_gene_depths(gene_table, depth_table):
    table = (
        gene_table.join(depth_table.set_index("seqid"), on="seqid")[
            lambda x: (x.pos >= x.left) & (x.pos <= x.right)
        ]
        .groupby(["seqid", "left", "right", "strand"])
        .depth.sum()
        .reset_index()
        .assign(
            length=lambda x: x.right - x.left, mean_depth=lambda x: x.depth / x.length
        )
        .rename(columns={"depth": "total_depth"})
    )
    return table


if __name__ == "__main__":
    gff_inpath = sys.argv[1]
    depth_inpath = sys.argv[2]
    outpath = sys.argv[3]

    gene_depth_tables = []
    with open(gff_inpath) as gff_stream, open(depth_inpath) as depth_stream:
        gff_chunk_iter = group_and_tabularize_genes_by_seqid(
            parse_gff_to_genes(tqdm(gff_stream))
        )
        depth_chunk_iter = group_and_tabularize_depths_by_seqid(
            parse_depths(depth_stream)
        )
        for seqid, gene_table, depth_table in join_chunk_iters(
            gff_chunk_iter, depth_chunk_iter
        ):
            gene_depth_tables.append(
                calculate_mean_gene_depths(gene_table, depth_table).assign(seqid=seqid)
            )
    breakpoint()
    gene_depth_tables = pd.concat(gene_depth_tables)

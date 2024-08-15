# {{{2 Data Configuration


config["figures"]["submission"] = []


config["genomes"] = pd.read_table(
    "meta/genomes.tsv", index_col=["species_id", "strain"]
)

config["mgen"] = pd.read_table("meta/mgen_to_reads.tsv", index_col="mgen_id")

config["mgen_group"] = (
    pd.read_table("meta/mgen_group.tsv")
    .groupby("mgen_group")
    .mgen_id.apply(list)
)

# {{{2 Data Configuration


config["figures"]["submission"] = []


config["genomes"] = pd.read_table(
    "meta/genomes.tsv", index_col=["species_id", "strain"]
)

config["mgen_group"] = (
    pd.read_table("meta/mgen_group.tsv", names=["mgen_id", "mgen_group"])
    .groupby("mgen_group")
    .mgen_id.apply(list)
)

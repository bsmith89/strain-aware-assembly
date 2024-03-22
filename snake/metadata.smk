# {{{2 Data Configuration


config["figures"]["submission"] = []


config["genomes"] = pd.read_table(
    "meta/genomes.tsv", index_col=["species_id", "strain"]
)

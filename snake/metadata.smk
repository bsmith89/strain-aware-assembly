# {{{2 Data Configuration


config["figures"]["submission"] = []


# config["genomes"] = pd.read_table(
#     "meta/genomes.tsv", index_col=["species_id", "strain"]
# )
#
config["mgen"] = pd.read_table("meta/mgen_to_reads.tsv", index_col="mgen_id")

config["mgen_group"] = (
    pd.read_table("meta/mgen_group.tsv").groupby("mgen_group").mgen_id.apply(list)
)

config["genome"] = pd.read_table("meta/genome.tsv", dtype=str, index_col=["genome_id"])

_genome_group = pd.read_table(
    "meta/genome_group.tsv",
    dtype=str,
)
for genome_group_id, d in _genome_group.groupby("genome_group_id"):
    config["genome_group"][genome_group_id] = list(d.genome_id)

config["simulated_community"] = pd.read_table("meta/simulated_community.tsv")

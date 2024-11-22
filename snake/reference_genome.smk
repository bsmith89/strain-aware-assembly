def genbank_genomic_ftp_url(accession, assembly):
    prefix = accession[:3]
    n1to3 = accession[4:7]
    n4to6 = accession[7:10]
    n7to9 = accession[10:13]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{n1to3}/{n4to6}/{n7to9}/{accession}_{assembly}/{accession}_{assembly}_genomic.fna.gz"


# e.g. https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/512/915/GCA_000512915.1_ASM51291v1/GCA_000512915.1_ASM51291v1_genomic.fna.gz


rule download_genbank_genome:
    output:
        "raw/genbank/{accession}_{assembly}.fn",
    params:
        url=lambda w: genbank_genomic_ftp_url(w.accession, w.assembly),
    shell:
        "curl '{params.url}' | zcat > {output}"


localrules:
    download_genbank_genome,


def alias_genbank_genome_input(w):
    accession, assembly = config["genomes"].loc[(w.species, w.strain)][
        ["genbank", "assembly"]
    ]
    return f"raw/genbank/{accession}_{assembly}.fn"


rule alias_genbank_genome:
    output:
        "data/genbank/{species}.{strain}.fn",
    input:
        alias_genbank_genome_input,
    shell:
        alias_recipe


rule standardize_contig_naming_project_reference_genome:
    output:
        "data/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        "sed '/>/s:>\(.*\):>{wildcards.genome}_\\1:' {input} > {output}"


rule normalize_genome_sequence:
    output:
        "{stem}.norm.fn",
    input:
        "{stem}.fn",
    shell:
        "sed '/^>/!s/[^ACGT]/N/g' {input} > {output}"


rule tile_reference_genome:
    output:
        "{stem}.tiles-k{ksize}.fn",
    input:
        script="scripts/tile_fasta.py",
        fn="{stem}.fn",
    wildcard_constraints:
        genome=noperiod_wc,
        ksize=integer_wc,
    params:
        length=lambda w: int(w.ksize),
        overlap=lambda w: int(w.ksize) - 1,
    shell:
        "{input.script} {params.length} {params.overlap} {input.fn} > {output}"


rule tile_reference_genome2:
    output:
        "{stem}.tiles-k{ksize}-o{overlap}.fn",
    input:
        script="scripts/tile_fasta.py",
        fn="{stem}.fn",
    wildcard_constraints:
        genome=noperiod_wc,
        ksize=integer_wc,
        overlap=integer_wc,
    params:
        length=lambda w: int(w.ksize),
        overlap=lambda w: int(w.overlap),
    shell:
        "{input.script} {params.length} {params.overlap} {input.fn} > {output}"


rule genome_fasta_to_fastq:
    """
    Convert a FASTA formatted file into FASTQ.
    \
    Input/output patterns are limited to files found in */genome/* in order to
    prevent circular dependencies.
    \
    """
    output:
        "{stemA}/genome/{stemB}.fq.gz",
    input:
        "{stemA}/genome/{stemB}.fn",
    conda:
        "conda/seqtk.yaml"
    shell:
        "seqtk seq -F '#' {input} | gzip -c > {output}"


rule alias_tiled_genome_as_reads:
    output:
        r1="data/reads/{genome}_tiles_k{ksize}/r1.tiles.fq.gz",
        r2="data/reads/{genome}_tiles_k{ksize}/r2.tiles.fq.gz",
    input:
        tiles="data/genome/{genome}.norm.tiles-k{ksize}.fq.gz",
    shell:
        """
        ln -rs {input.tiles} {output.r1}
        ln -rs {input.tiles} {output.r2}
        """


ruleorder: alias_tiled_genome_as_reads > alias_cleaned_reads


rule simulate_wgs_reads:
    output:
        r1="data/reads/{genome}_sim_len150_seed{seed}_cov{cov}/r1.sim.fq.gz",
        r2="data/reads/{genome}_sim_len150_seed{seed}_cov{cov}/r2.sim.fq.gz",
    input:
        fasta="data/genome/{genome}.fn",
    params:
        coverage=lambda w: int(w.cov) / 100,
        outdir="data/reads/{genome}_sim_len150_seed{seed}_cov{cov}",
    conda:
        "conda/art_read_sim.yaml"
    shell:
        """
        art_illumina --rndSeed {wildcards.seed} \
                --seqSys HS25 \
                --len 150 --paired --mflen 450 --sdev 30 \
                --fcov {params.coverage} \
                --noALN \
                --in {input.fasta} --out {params.outdir}/
        gzip -c {params.outdir}/1.fq > {output.r1} && rm {params.outdir}/1.fq
        gzip -c {params.outdir}/2.fq > {output.r2} && rm {params.outdir}/2.fq
        """


rule combine_community_wgs_reads:
    output:
        r1="data/reads/{community}_simcom_len150_seed{seed}_cov{cov}/r1.sim.fq.gz",
        r2="data/reads/{community}_simcom_len150_seed{seed}_cov{cov}/r2.sim.fq.gz",
    input:
        r1=lambda w: [
            "data/reads/{genome}_sim_len150_seed{seed}_cov{cov}/r1.sim.fq.gz".format(
                genome=sim_genome.genome_id,
                seed=w.seed,
                cov=int(int(w.cov) * sim_genome.fraction),
            )
            for _, sim_genome in config["simulated_community"][
                lambda x: x.community_id == w.community
            ].iterrows()
        ],
        r2=lambda w: [
            "data/reads/{genome}_sim_len150_seed{seed}_cov{cov}/r2.sim.fq.gz".format(
                genome=sim_genome.genome_id,
                seed=w.seed,
                cov=int(int(w.cov) * sim_genome.fraction),
            )
            for _, sim_genome in config["simulated_community"][
                lambda x: x.community_id == w.community
            ].iterrows()
        ],
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """


rule construct_mgen_group_input_table:
    output:
        "data/group/{group}/r.proc.kmtricks_input.txt",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.proc.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.proc.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        entry_pattern="$mgen : data/reads/$mgen/r1.proc.fq.gz ; data/reads/$mgen/r2.proc.fq.gz",
        mgen_list=lambda w: config["mgen_group"][w.group],
    shell:
        """
        for mgen in {params.mgen_list}; do echo "{params.entry_pattern}"; done > {output}
        """


rule construct_three_genome_input_table:
    output:
        "data/three_genomes.kmtricks_input.txt",
    input:
        ecoli_mg1655="data/genbank/ecoli.mg1655.fn",
        ecoli_o121h19="data/genbank/ecoli.o121h19.fn",
        bdorei_dsm17855="data/genbank/bdorei.dsm17855.fn",
    shell:
        dd(
            """
        cat <<EOF > {output}
        ecoli_mg1655 : data/genbank/ecoli.mg1655.fn
        ecoli_o121h19 : data/genbank/ecoli.o121h19.fn
        bdorei_dsm17855 : data/genbank/bdorei.dsm17855.fn
        EOF
        """
        )


rule construct_two_genome_input_table:
    output:
        "data/two_genomes.kmtricks_input.txt",
    input:
        ecoli_mg1655="data/genbank/ecoli.mg1655.fn",
        ecoli_o121h19="data/genbank/ecoli.o121h19.fn",
    shell:
        dd(
            """
        cat <<EOF > {output}
        ecoli_mg1655 : data/genbank/ecoli.mg1655.fn
        ecoli_o121h19 : data/genbank/ecoli.o121h19.fn
        EOF
        """
        )

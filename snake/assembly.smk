rule start_shell_ggcat:
    container:
        config["container"]["ggcat"]
    shell:
        "bash"


rule start_conda_shell:
    output:
        "start_shell.{conda}",
    conda:
        "conda/{conda}.yaml"
    shell:
        "bash"


rule start_conda_ipython:
    output:
        "start_ipython.{conda}",
    conda:
        "conda/{conda}.yaml"
    shell:
        "ipython"


#
#
# use rule start_shell as start_shell_ggcat with:
#     container:
#         config["container"]["ggcat"]
#
#
# use rule start_shell as start_shell_kmtricks with:
#     conda:
#         "conda/kmtricks.yaml"
#


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
                cov=int(int(w.cov) * sim_genome.fraction),  # NOTE (2024-11-19): Metadata has coverage in correct units. Filenames for this and single genomes are 100x.
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


rule run_kmtricks_repart:
    output:
        repart_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/repart.flag",
    input:
        "{stem}.kmtricks_input.txt",
    params:
        workdir="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
        ksize=lambda w: int(w.ksize),
        num_partitions=32,  # NOTE: This determines the parallelizability of the merging step.
    conda:
        "conda/kmtricks.yaml"
    threads: 8
    shell:
        """
        rm -r {params.workdir}  # Snakemake auto-builds this directory, so I need to delete it for kmtricks to run.
        kmtricks repart --run-dir {params.workdir} \
                --threads {threads} --verbose info \
                --kmer-size {params.ksize} \
                --nb-partitions {params.num_partitions} \
                --file {input}
        mkdir {params.workdir}/checkpoints
        touch {output.repart_flag}
        """


rule run_kmtricks_superk:
    output:
        superk_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/{mgen}.superk.flag",
    input:
        repart_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/repart.flag",
    params:
        workdir="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
    conda:
        "conda/kmtricks.yaml"
    threads: 1
    shell:
        """
        kmtricks superk --threads {threads} --verbose info --run-dir {params.workdir} --id {wildcards.mgen}
        touch {output.superk_flag}
        """


rule run_kmtricks_count:
    output:
        count_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/{mgen}.count.flag",
    input:
        superk_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/{mgen}.superk.flag",
    params:
        workdir="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
    conda:
        "conda/kmtricks.yaml"
    threads: 8
    shell:
        """
        kmtricks count --run-dir {params.workdir} \
                --threads {threads} --verbose info \
                --mode kmer \
                --hard-min 0 \
                --id {wildcards.mgen}
        touch {output.count_flag}
        """


# NOTE: --hard-min 0 and --share-min 1 always
# or else some 0s will be censoring instead of actually
# missing. (That's because counts below --hard-min are discarded
# from individual samples and counts below --soft-min can be
# discarded from individual samples if they're solid in fewer than
# --share-min samples.
# With --hard-min 0 and --share-min 1, counts for a kmer are either
# kept entirely or discarded entirely.
rule run_kmtricks_merge:
    output:
        merge_flag="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/matrices/matrix_{part}.count.txt",
    input:
        count_flags=lambda w: [
            f"data/group/{w.group}/{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/checkpoints/{mgen}.count.flag"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        workdir="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
        mincount=lambda w: int(w.mincount),
        recurrence=lambda w: int(w.recurrence),
    conda:
        "conda/kmtricks.yaml"
    threads: 1
    shell:
        """
        kmtricks merge --run-dir {params.workdir} \
            --threads 1 --verbose info \
            --mode kmer:count:text \
            --share-min 1 \
            --soft-min {params.mincount} \
            --recurrence-min {params.recurrence} \
            --partition-id {wildcards.part}
        touch {output.merge_flag}
        """


rule drop_short_tips_from_ggcat_fasta:
    output:
        "{stem}.ggcat-k{ksize}-{unitig_source}.droptips-{length}.fn",
    wildcard_constraints:
        length=single_param_wc,
    input:
        script="scripts/drop_short_tips_from_ggcat_fasta.py",
        fasta="{stem}.ggcat-k{ksize}-{unitig_source}.fn",
    params:
        length_mult=lambda w: float(w.length),
        ksize=lambda w: int(w.ksize),
    conda:
        "conda/strainzip.yaml"
    shell:
        "{input.script} {params.ksize} {params.length_mult} {input.fasta} {output}"


rule rerun_ggcat_on_droptips_fasta:
    output:
        "{stem}.ggcat-k{ksize}-{unitig_source}-droptips.fn",
    input:
        "{stem}.ggcat-k{ksize}-{unitig_source}.droptips-2.fn",
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d) && echo $input_dir
        ln -s $(realpath {input}) $input_dir/$(basename {input}).fa
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_kmc_on_reads:
    output:
        pre="{stemA}/r.{stemB}.kmc-k{ksize}.kmc_pre",
        suf="{stemA}/r.{stemB}.kmc-k{ksize}.kmc_suf",
    wildcard_constraints:
        ksize=integer_wc,
    input:
        r1="{stemA}/r1.{stemB}.fq.gz",
        r2="{stemA}/r2.{stemB}.fq.gz",
    params:
        ksize=lambda w: int(w.ksize),
        dbname="{stemA}/r.{stemB}.kmc-k{ksize}",
    threads: 1
    conda:
        "conda/kmc.yaml"
    shell:
        """
        workdir=$(mktemp -d)

        # Construct input file list.
        filelist=$(mktemp)
        for file in {input}
        do
            echo $file
        done > $filelist

        echo $workdir $filelist

        kmc -v -k{params.ksize} -fq -ci1 -cs1000000000 -t{threads} @$filelist {params.dbname} $workdir

        # Cleanup
        rm -r $filelist $workdir
        """


rule build_kmc_mask_from_fasta:
    output:
        pre="{stem}.kmc-k{ksize}-mask.kmc_pre",
        suf="{stem}.kmc-k{ksize}-mask.kmc_suf",
    input:
        fasta="{stem}.fn",
    params:
        ksize=lambda w: int(w.ksize),
        dbname="{stem}.kmc-k{ksize}-mask",
    threads: 8
    conda:
        "conda/kmc.yaml"
    shell:
        """
        workdir=$(mktemp -d)
        echo $workdir

        kmc -v -k{params.ksize} -fa -ci1 -cs2 -t{threads} {input.fasta} {params.dbname} $workdir

        # Cleanup
        rm -r $workdir
        """


rule alias_kmc_counts_denovo0:
    """
    NOTE (2024-11-15): This rule is only necessary because the file naming
    changed so much when using the "droptips" masking approach below.
    """
    output:
        "data/group/{group}/reads/{mgen}/r.{stem}.kmc-k{ksize}-denovo0.kmc_{pre_or_suf}",
    input:
        reads_pre="data/reads/{mgen}/r.{stem}.kmc-k{ksize}.kmc_{pre_or_suf}",
    shell:
        alias_recipe


rule filter_kmc_counts_by_notips_contigs:
    output:
        pre="data/group/{group}/reads/{mgen}/r.{stem}.kmc-k{ksize}-{unitig_source}-droptips.kmc_pre",
        suf="data/group/{group}/reads/{mgen}/r.{stem}.kmc-k{ksize}-{unitig_source}-droptips.kmc_suf",
    input:
        reads_pre="data/reads/{mgen}/r.{stem}.kmc-k{ksize}.kmc_pre",
        reads_suf="data/reads/{mgen}/r.{stem}.kmc-k{ksize}.kmc_suf",
        mask_pre="data/group/{group}/r.{stem}.ggcat-k{ksize}-{unitig_source}-droptips.kmc-k{ksize}-mask.kmc_pre",
        mask_suf="data/group/{group}/r.{stem}.ggcat-k{ksize}-{unitig_source}-droptips.kmc-k{ksize}-mask.kmc_suf",
    params:
        reads_dbname="data/reads/{mgen}/r.{stem}.kmc-k{ksize}",
        mask_dbname="data/group/{group}/r.{stem}.ggcat-k{ksize}-{unitig_source}-droptips.kmc-k{ksize}-mask",
        out_dbname="data/group/{group}/reads/{mgen}/r.{stem}.kmc-k{ksize}-{unitig_source}-droptips",
        ksize=lambda w: int(w.ksize),
    conda:
        "conda/kmc.yaml"
    shell:
        """
        kmc_tools simple {params.mask_dbname} {params.reads_dbname} intersect {params.out_dbname} -ocright
        """


# NOTE: This output is very large and this rule may be only
# useful for debugging. In reality, I should only stream the output.
# e.g. with `<(kmc_tools transform  data/group/{w.group}/reads/{mgen}/r.{w.stem} dump >(cat))`
rule dump_kmc_counts:
    output:
        "{stem}.kcounts.tsv",
    input:
        pre="{stem}.kmc_pre",
        suf="{stem}.kmc_suf",
    params:
        dbname="{stem}",
    conda:
        "conda/kmc.yaml"
    shell:
        "kmc_tools transform {params.dbname} dump -s {output}"


rule merge_kmc_counts_to_table:
    output:
        "data/group/{group}/r.{stem}.kcounts_merged.tsv",
    input:
        pre_and_suf=lambda w: [
            f"data/group/{w.group}/reads/{mgen}/r.{w.stem}.kmc_{pre_or_suf}"
            for mgen, pre_or_suf in product(
                config["mgen_group"][w.group], ["pre", "suf"]
            )
        ],
    params:
        header=lambda w: "kmer\t" + "\t".join(config["mgen_group"][w.group]),
        args=lambda w: [
            f"<(kmc_tools transform  data/group/{w.group}/reads/{mgen}/r.{w.stem} dump >(cat))"
            for mgen in config["mgen_group"][w.group]
        ],
    conda:
        "conda/kmc.yaml"
    shell:
        """
        echo "{params.header}" > {output}
        scripts/merge_kmc_counts.py {params.args} | tqdm >> {output}
        """


rule load_kmc_merged_table_to_sqlite:
    output:
        "{stem}.kmc-k{ksize}-{unitig_source}-droptips.db",
    input:
        script="scripts/sample_list_to_sqlite_table_script.py",
        counts="{stem}.kmc-k{ksize}-{unitig_source}-droptips.kcounts_merged.tsv",
        sample_list="{stem}.kmtricks_input.txt",
    shell:
        """
        tmpdb=$(mktemp) && echo $tmpdb
        {input.script} {input.sample_list} {wildcards.ksize} | sqlite3 $tmpdb
        cat {input.counts} | sed '1,1d' | tqdm --unit-scale 1 | sqlite3 -separator '\t' $tmpdb '.import /dev/stdin count_'
        mv $tmpdb {output}
        """


rule load_kmtricks_output_to_sqlite:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    input:
        script="scripts/sample_list_to_sqlite_table_script.py",
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
        sample_list="{stem}.kmtricks_input.txt",
    shell:
        """
        tmpdb=$(mktemp) && echo $tmpdb
        {input.script} {input.sample_list} {wildcards.ksize} | sqlite3 $tmpdb
        cat {input.counts} | tqdm --unit-scale 1 | sqlite3 -separator ' ' $tmpdb '.import /dev/stdin count_'
        mv $tmpdb {output}
        """


rule run_ggcat_on_kmtricks_kmers:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-denovo.fn",
    input:
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        for file in {input.counts}
        do
            fasta=$input_dir/$(basename $file).fifo.fa
            echo making fifo $fasta
            mkfifo $fasta
            awk -v OFS="\\n" '{{print ">"NR,$1}}' $file > $fasta &
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_ggcat_on_kmtricks_kmers_and_include_xjin_refs:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-withref.fn",
    input:
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"]["xjin"]
        ],
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        for file in {input.refs}
        do
            ln -s $(realpath $file) $input_dir/$(basename $file).fa
        done

        for file in {input.counts}
        do
            fasta=$input_dir/$(basename $file).fifo.fa
            echo making fifo $fasta
            mkfifo $fasta
            awk -v OFS="\\n" '{{print ">"NR,$1}}' $file > $fasta &
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_ggcat_on_kmtricks_kmers_and_include_megahit_contigs:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-withmegahit.fn",
    input:
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
        contigs="{stem}.megahit-full-k{ksize}-unfiltered.fn",
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        # MEGAHIT contigs
        ln -rs {input.contigs} $input_dir/$(basename {input.contigs}).fa

        # Kmtricks counts
        for file in {input.counts}
        do
            fasta=$input_dir/$(basename $file).fifo.fa
            echo making fifo $fasta
            mkfifo $fasta
            awk -v OFS="\\n" '{{print ">"NR,$1}}' $file > $fasta &
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_ggcat_on_xjin_refs:
    output:
        "data/xjin_ref.ggcat-k{ksize}-refonly.fn",
    input:
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"]["xjin"]
        ],
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        for file in {input.refs}
        do
            ln -s $(realpath $file) $input_dir/$(basename $file).fa
        done

        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_ggcat_on_reads:
    output:
        "data/group/{group}/r.{stem}.ggcat-k{ksize}-denovo.fn",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        ggcat build -s 2 -j {threads} -o {output} -k {wildcards.ksize} -e {input.r1} {input.r2}
        """


rule run_ggcat_on_reads_no_min:
    output:
        "data/group/{group}/r.{stem}.ggcat-k{ksize}-denovo0.fn",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e {input.r1} {input.r2}
        """


rule run_ggcat_on_reads_and_include_megahit_contigs:
    output:
        "data/group/{group}/r.{stem}.ggcat-k{ksize}-withmegahit.fn",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        contigs="data/group/{group}/r.{stem}.megahit-full-k{ksize}-unfiltered.fn",
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        # MEGAHIT contigs
        ln -rs {input.contigs} $input_dir/megahit_contigs.fa

        ggcat build -s 2 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/megahit_contigs.fa {input.r1} {input.r2}

        rm -r $input_dir
        """


rule run_ggcat_on_reads_and_include_megahit_contigs_twice:
    output:
        "data/group/{group}/r.{stem}.ggcat-k{ksize}-withmegahit2.fn",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        contigs="data/group/{group}/r.{stem}.megahit-full-k{ksize}-unfiltered.fn",
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        # MEGAHIT contigs
        ln -rs {input.contigs} $input_dir/megahit_contigs1.fa
        ln -rs {input.contigs} $input_dir/megahit_contigs2.fa

        ggcat build -s 2 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/* {input.r1} {input.r2}

        rm -r $input_dir
        """


rule run_ggcat_on_reads_and_include_megahit_contigs_min3:
    output:
        "data/group/{group}/r.{stem}.ggcat-k{ksize}-withmegahit3.fn",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        contigs="data/group/{group}/r.{stem}.megahit-full-k{ksize}-unfiltered.fn",
    container:
        config["container"]["ggcat"]
    threads: 48
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        # MEGAHIT contigs
        ln -rs {input.contigs} $input_dir/megahit_contigs1.fa
        ln -rs {input.contigs} $input_dir/megahit_contigs2.fa
        ln -rs {input.contigs} $input_dir/megahit_contigs3.fa

        ggcat build -s 3 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/* {input.r1} {input.r2}

        rm -r $input_dir
        """


rule calculate_mean_unitig_depths_across_samples_from_kmtricks:
    output:
        "data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.unitig_depth.nc",
    wildcard_constraints:
        unitig_source=noperiod_wc,
    input:
        fasta="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.fn",
        db="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    conda:
        "conda/strainzip.yaml"
    params:
        sample_list=lambda w: ",".join(config["mgen_group"][w.group]),
    threads: 48
    shell:
        """
        tmpdb=$(mktemp) && echo $tmpdb
        strainzip depth --verbose --preload --tmpdb $tmpdb -p {threads} {input.fasta} {input.db} {wildcards.ksize} {params.sample_list} {output}
        rm $tmpdb
        """


# rule calculate_mean_unitig_depths_across_samples_from_kmc_droptips:
#     output:
#         "data/group/{group}/{stem}.ggcat-k{ksize}-{unitig_source}-droptips.unitig_depth.nc",
#     wildcard_constraints:
#         unitig_source=noperiod_wc,
#     input:
#         fasta="data/group/{group}/{stem}.ggcat-k{ksize}-{unitig_source}-droptips.fn",
#         db="data/group/{group}/{stem}.kmc-k{ksize}-{unitig_source}-droptips.db",
#     conda:
#         "conda/strainzip.yaml"
#     params:
#         sample_list=lambda w: ",".join(config["mgen_group"][w.group]),
#     threads: 48
#     shell:
#         """
#         tmpdb=$(mktemp) && echo $tmpdb
#         strainzip depth --verbose --preload --tmpdb $tmpdb -p {threads} {input.fasta} {input.db} {wildcards.ksize} {params.sample_list} {output}
#         rm $tmpdb
#         """


rule calculate_mean_unitig_depths_across_samples_from_kmc_one_step:
    output:
        "data/group/{group}/{stem}.ggcat-k{ksize}-{unitig_source}.unitig_depth.nc",
    wildcard_constraints:
        ksize=integer_wc,
    input:
        merge_script="scripts/merge_kmc_counts.py",
        mean_script="scripts/mean_unitig_kmer_depth.py",
        load_netcdf_script="scripts/load_unitig_depths_to_netcdf.py",
        fasta="data/group/{group}/{stem}.ggcat-k{ksize}-{unitig_source}.fn",
        pre_and_suf=lambda w: [
            f"data/group/{w.group}/reads/{mgen}/{w.stem}.kmc-k{w.ksize}-{w.unitig_source}.kmc_{pre_or_suf}"
            for mgen, pre_or_suf in product(
                config["mgen_group"][w.group], ["pre", "suf"]
            )
        ],
    params:
        ksize=lambda w: int(w.ksize),
        sample_names=lambda w: ",".join(config["mgen_group"][w.group]),
        args=lambda w: [
            f"<(scripts/kmc_dump_to_stdout.sh data/group/{w.group}/reads/{mgen}/{w.stem}.kmc-k{w.ksize}-{w.unitig_source})"
            for mgen in config["mgen_group"][w.group]
        ],
    conda:
        "conda/strainzip_kmc.yaml"
    shell:
        """
        {input.merge_script} {params.args} \
                | {input.mean_script} {params.ksize} {params.sample_names} {input.fasta} {output}
        """


# rule load_ggcat_without_depths_to_sz:
#     output:
#         "{stem}.ggcat-k{ksize}-{unitig_source}.nodepth.sz",
#     wildcard_constraints:
#         unitig_source=single_param_wc,
#     input:
#         fasta="{stem}.ggcat-k{ksize}-{unitig_source}.fn",
#     conda:
#         "conda/strainzip.yaml"
#     shell:
#         """
#         strainzip load --verbose {wildcards.ksize} {input.fasta} {output}
#         """


rule load_ggcat_with_kmc_depths_to_sz:
    output:
        "{stem}.ggcat-k{ksize}-{unitig_source}.sz",
    wildcard_constraints:
        unitig_source=noperiod_wc,
        ksize=integer_wc,
    input:
        fasta="{stem}.ggcat-k{ksize}-{unitig_source}.fn",
        depth="{stem}.ggcat-k{ksize}-{unitig_source}.unitig_depth.nc",
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip load --verbose {wildcards.ksize} --depth {input.depth} {input.fasta} {output}
        """


# # Because the masking only affects the tips.
# rule rename_tip_trimmed_kmc_graph_to_drop_masking_infix:
#     output:
#         "{stem}.ggcat-k{ksize}-{unitig_source}.notips-2.sz"
#     input:
#         "{stem}.ggcat-k{ksize}-{unitig_source}-droptips.notips-2.sz"
#     shell:
#         alias_recipe
#
#
# ruleorder: rename_tip_trimmed_kmc_graph_to_drop_masking_infix > trim_tips


rule load_ggcat_with_kmtricks_depths_to_sz:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.sz",
    wildcard_constraints:
        unitig_source=single_param_wc,
    input:
        fasta="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.fn",
        depth="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.unitig_depth.nc",
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip load --verbose {wildcards.ksize} --depth {input.depth} {input.fasta} {output}
        """


rule trim_tips:
    output:
        "{stem}.notips-{length}.sz",
    wildcard_constraints:
        length=single_param_wc,
    input:
        "{stem}.sz",
    params:
        length=lambda w: float(w.length),
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip trim --verbose -p {threads} --num-kmer-lengths {params.length} {input} {output}"


# FIXME (2024-05-20): Drop this when I've updated the CLI permanently
rule trim_tips_old:
    output:
        "{stem}.notips.sz",
    wildcard_constraints:
        length=single_param_wc,
    input:
        "{stem}.sz",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip trim --verbose {input} {output}"


rule trim_tips_unpressed:
    output:
        "{stem}.notips-{length}-unpressed.sz",
    input:
        "{stem}.sz",
    params:
        length=lambda w: float(w.length),
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip trim --verbose -p {threads} --no-press --num-kmer-lengths {params.length} {input} {output}"


rule smooth_depths:
    output:
        "{stem}.smoothed-{eps}.sz",
    wildcard_constraints:
        eps=integer_wc,
    input:
        "{stem}.sz",
    params:
        eps=lambda w: 10 ** (-int(w.eps)),
    conda:
        "conda/strainzip.yaml"
    threads: 48
    shell:
        """
        strainzip smooth --verbose -p {threads} --eps {params.eps} {input} {output}
        """


rule smooth_depths_enforce_symmetry:
    output:
        "{stem}.smoothed-sym-{eps}.sz",
    wildcard_constraints:
        eps=integer_wc,
    input:
        "{stem}.sz",
    params:
        eps=lambda w: 10 ** (-int(w.eps)),
    conda:
        "conda/strainzip.yaml"
    threads: 48
    shell:
        """
        strainzip smooth --verbose -p {threads} --eps {params.eps} --enforce-symmetry {input} {output}
        """


rule unzip_safe_only_junctions:
    output:
        final="{stem}.unzip_safe-{model}-{thresh}-{rounds}.sz",
    wildcard_constraints:
        model=single_param_wc,
        thresh=single_param_wc,
        rounds=single_param_wc,
    input:
        "{stem}.sz",
    log:
        checkpoint_dir=directory(
            "{stem}.unzip_safe-{model}-{thresh}-{rounds}.checkpoints.d"
        ),
    params:
        model=lambda w: {
            "lognorm2": "OffsetLogNormal",
            "norm": "Normal --model-hyperparameters tol=1e-4",
            "normscaled": "NormalScaled --model-hyperparameters alpha=0.5",
            "lapl": "Laplace",
            "t5": "StudentsT --model-hyperparameters df=5",
            "huber": "Huber --model-hyperparameters delta=1",
        }[w.model],
        min_depth=1.0,
        score_thresh=lambda w: float(w.thresh),
        relative_error_thresh=0.1,
        absolute_error_thresh=1.0,
        max_rounds=lambda w: int(w.rounds),
        excess_thresh=0,
        completeness_thresh=1,
    conda:
        "conda/strainzip.yaml"
    threads: 48
    shell:
        """
        # FIXME: Figure out why setting environmental variables here is necessary.
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1 --xla_force_host_platform_device_count=8"
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREAD=1
        export NUM_INTER_THREADS=1
        export NUM_INTRA_THREADS=1

        mkdir -p {log.checkpoint_dir}

        strainzip unzip --verbose -p {threads} \
                --min-depth {params.min_depth} \
                --skip-canonical --skip-large --skip-extra-large \
                --max-rounds {params.max_rounds} --model {params.model} \
                --score aic --score-thresh {params.score_thresh} \
                --relative-error-thresh {params.relative_error_thresh} \
                --absolute-error-thresh {params.absolute_error_thresh} \
                --excess-thresh {params.excess_thresh} \
                --completeness-thresh {params.completeness_thresh} \
                --checkpoint-dir {log.checkpoint_dir} \
                {input} {output.final}
        """


rule unzip_junctions:
    output:
        final="{stem}.unzip-{model}-{thresh}-{rounds}.sz",
    wildcard_constraints:
        model=single_param_wc,
        thresh=single_param_wc,
        rounds=single_param_wc,
    input:
        "{stem}.sz",
    log:
        checkpoint_dir=directory("{stem}.unzip-{model}-{thresh}-{rounds}.checkpoints.d"),
    params:
        model=lambda w: {
            "lognorm2": "OffsetLogNormal",
            "norm": "Normal --model-hyperparameters tol=1e-4",
            "normscaled": "NormalScaled --model-hyperparameters alpha=0.5",
            "lapl": "Laplace",
            "t5": "StudentsT --model-hyperparameters df=5",
            "huber": "Huber --model-hyperparameters delta=1",
        }[w.model],
        min_depth=1.0,
        score_thresh=lambda w: float(w.thresh),
        relative_error_thresh=0.1,
        absolute_error_thresh=1.0,
        max_rounds=lambda w: int(w.rounds),
        excess_thresh=0,
        completeness_thresh=1,
    conda:
        "conda/strainzip.yaml"
    threads: 48
    shell:
        """
        # FIXME: Figure out why setting environmental variables here is necessary.
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1 --xla_force_host_platform_device_count=8"
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREAD=1
        export NUM_INTER_THREADS=1
        export NUM_INTRA_THREADS=1

        mkdir -p {log.checkpoint_dir}

        strainzip unzip --verbose -p {threads} \
                --min-depth {params.min_depth} \
                --skip-extra-large --max-rounds {params.max_rounds} --model {params.model} \
                --score aic --score-thresh {params.score_thresh} \
                --relative-error-thresh {params.relative_error_thresh} \
                --absolute-error-thresh {params.absolute_error_thresh} \
                --excess-thresh {params.excess_thresh} \
                --completeness-thresh {params.completeness_thresh} \
                --checkpoint-dir {log.checkpoint_dir} \
                {input} {output.final}
        """


rule unzip_junctions_no_balancing:
    output:
        final="{stem}.unzip-{model}-nobal-{thresh}-{rounds}.sz",
    wildcard_constraints:
        model=single_param_wc,
        thresh=single_param_wc,
        rounds=single_param_wc,
    input:
        "{stem}.sz",
    log:
        checkpoint_dir=directory(
            "{stem}.unzip-{model}-nobal-{thresh}-{rounds}.checkpoints.d"
        ),
    params:
        model=lambda w: {
            "lognorm2": "OffsetLogNormal",
            "norm": "Normal --model-hyperparameters tol=1e-4",
            "normscaled": "NormalScaled --model-hyperparameters alpha=0.5",
            "lapl": "Laplace",
            "t5": "StudentsT --model-hyperparameters df=5",
            "huber": "Huber --model-hyperparameters delta=1",
        }[w.model],
        min_depth=1.0,
        score_thresh=lambda w: float(w.thresh),
        relative_error_thresh=0.1,
        absolute_error_thresh=1.0,
        max_rounds=lambda w: int(w.rounds),
        excess_thresh=0,
        completeness_thresh=1,
    conda:
        "conda/strainzip.yaml"
    threads: 48
    shell:
        """
        # FIXME: Figure out why setting environmental variables here is necessary.
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1 --xla_force_host_platform_device_count=8"
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREAD=1
        export NUM_INTER_THREADS=1
        export NUM_INTRA_THREADS=1

        mkdir -p {log.checkpoint_dir}

        strainzip unzip --verbose -p {threads} \
                --min-depth {params.min_depth} \
                --skip-extra-large --max-rounds {params.max_rounds} --model {params.model} \
                --score aic --score-thresh {params.score_thresh} \
                --no-balance \
                --relative-error-thresh {params.relative_error_thresh} \
                --absolute-error-thresh {params.absolute_error_thresh} \
                --excess-thresh {params.excess_thresh} \
                --completeness-thresh {params.completeness_thresh} \
                --checkpoint-dir {log.checkpoint_dir} \
                {input} {output.final}
        """


rule unzip_junctions_no_balancing_no_low_depth:
    output:
        final="{stem}.unzip-{model}-nobal-nocull-{thresh}-{rounds}.sz",
    wildcard_constraints:
        model=single_param_wc,
        thresh=single_param_wc,
        rounds=single_param_wc,
    input:
        "{stem}.sz",
    log:
        checkpoint_dir=directory(
            "{stem}.unzip-{model}-nobal-nocull-{thresh}-{rounds}.checkpoints.d"
        ),
    params:
        model=lambda w: {
            "lognorm2": "OffsetLogNormal",
            "norm": "Normal --model-hyperparameters tol=1e-4",
            "normscaled": "NormalScaled --model-hyperparameters alpha=0.5",
            "lapl": "Laplace",
            "t5": "StudentsT --model-hyperparameters df=5",
            "huber": "Huber --model-hyperparameters delta=1",
        }[w.model],
        score_thresh=lambda w: float(w.thresh),
        relative_error_thresh=0.1,
        absolute_error_thresh=1.0,
        max_rounds=lambda w: int(w.rounds),
        excess_thresh=0,
        completeness_thresh=1,
    conda:
        "conda/strainzip.yaml"
    threads: 48
    shell:
        """
        # FIXME: Figure out why setting environmental variables here is necessary.
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1 --xla_force_host_platform_device_count=8"
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREAD=1
        export NUM_INTER_THREADS=1
        export NUM_INTRA_THREADS=1

        mkdir -p {log.checkpoint_dir}

        strainzip unzip --verbose -p {threads} \
                --no-drop-low-depth \
                --skip-extra-large --max-rounds {params.max_rounds} --model {params.model} \
                --score aic --score-thresh {params.score_thresh} \
                --no-balance \
                --relative-error-thresh {params.relative_error_thresh} \
                --absolute-error-thresh {params.absolute_error_thresh} \
                --excess-thresh {params.excess_thresh} \
                --completeness-thresh {params.completeness_thresh} \
                --checkpoint-dir {log.checkpoint_dir} \
                {input} {output.final}
        """


rule benchmark_depth_model:
    output:
        "{stem}.model_benchmark-{model}-{thresh}.tsv",
    wildcard_constraints:
        model=single_param_wc,
        thresh=single_param_wc,
    input:
        "{stem}.pkl",
    params:
        model=lambda w: {
            "lognorm2": "OffsetLogNormal",
            "norm": "Normal --model-hyperparameters tol=1e-4",
            "normscaled": "NormalScaled --model-hyperparameters alpha=0.5",
            "lapl": "Laplace",
            "t5": "StudentsT --model-hyperparameters df=5",
            "huber": "Huber --model-hyperparameters delta=1",
        }[w.model],
        score_thresh=lambda w: float(w.thresh),
        relative_error_thresh=0.1,
        absolute_error_thresh=1.0,
        excess_thresh=0,
        completeness_thresh=1,
    conda:
        "conda/strainzip.yaml"
    threads: 36
    shell:
        """
        # FIXME: Figure out why setting environmental variables here is necessary.
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1 --xla_force_host_platform_device_count=8"
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREAD=1
        export NUM_INTER_THREADS=1
        export NUM_INTRA_THREADS=1

        strainzip benchmark --verbose -p {threads} \
                --model {params.model} \
                --score aic --score-thresh {params.score_thresh} \
                --relative-error-thresh {params.relative_error_thresh} \
                --absolute-error-thresh {params.absolute_error_thresh} \
                --excess-thresh {params.excess_thresh} \
                --completeness-thresh {params.completeness_thresh} \
                {input} {output}
        """


rule precluster_vertices:
    output:
        vertex="{stem}.preclust-e{exponent}-n{num_preclust}.vertex.tsv",
    input:
        "{stem}.sz",
    params:
        exponent=lambda w: int(w.exponent) / 100,
        num_preclust=lambda w: int(w.num_preclust),
    threads: 12
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip precluster --verbose {input} --exponent {params.exponent} {params.num_preclust} {output}
        """


rule cluster_vertices:
    output:
        vertex="{stem}.clust-e{exponent}-n{num}-d{thresh}.vertex.tsv",
        segment="{stem}.clust-e{exponent}-n{num}-d{thresh}.segment.tsv",
        depth="{stem}.clust-e{exponent}-n{num}-d{thresh}.depth.tsv",
        shared="{stem}.clust-e{exponent}-n{num}-d{thresh}.shared.tsv",
        meta="{stem}.clust-e{exponent}-n{num}-d{thresh}.meta.tsv",
    input:
        graph="{stem}.sz",
        preclust="{stem}.preclust-e{exponent}-n{num}.vertex.tsv",
    params:
        thresh=lambda w: int(w.thresh) / 1000,
        exponent=lambda w: int(w.exponent) / 100,
    threads: 12
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip cluster --verbose \
                --exponent {params.exponent} \
                {input.graph} {input.preclust} {params.thresh} {output.vertex} {output.segment} {output.depth} {output.shared} {output.meta}
        """


rule extract_unassembled_cluster_subgraph:
    output:
        graph="{stemA}.ggcat-{unitig_source}.{stemB}.clust-e{exponent}-d{thresh}.unassembled-c{clust}-r{radius}.sz",
    wildcard_constraints:
        radius=integer_wc,
        # unitig_source=single_param_wc,
    input:
        graph="{stemA}.ggcat-{unitig_source}.sz",
        segment="{stemA}.ggcat-{unitig_source}.{stemB}.clust-e{exponent}-d{thresh}.segment.tsv",
    params:
        clust=lambda w: int(w.clust),
        radius=lambda w: int(w.radius),
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip focus --verbose \
                --radius {params.radius} \
                --segments <(awk '$2 == {params.clust} {{print $1}}' {input.segment}) \
                {input.graph} \
                {output.graph}
        """


rule extract_assembled_cluster_subgraph:
    output:
        graph="{stem}.clust-{thresh}.assembled-c{clust}-r{radius}.sz",
    wildcard_constraints:
        radius=integer_wc,
    input:
        graph="{stem}.sz",
        vertex="{stem}.clust-{thresh}.vertex.tsv",
    params:
        clust=lambda w: int(w.clust),
        radius=lambda w: int(w.radius),
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip focus --verbose \
                --radius {params.radius} \
                --vertices <(awk '$2 == {params.clust} {{print $1}}' {input.vertex}) \
                {input.graph} \
                {output.graph}
        """


rule draw_graph:
    output:
        "{stem}.graph.pdf",
    input:
        "{stem}.sz",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip draw {input} {output}"


rule dump_assembly_contigs:
    output:
        fasta="{stemA}.ggcat-{unitig_source}.{stemB}.fn",
    wildcard_constraints:
        unitig_source=noperiod_wc,
    input:
        graph="{stemA}.ggcat-{unitig_source}.{stemB}.sz",
        fasta="{stemA}.ggcat-{unitig_source}.fn",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip dump_contigs --verbose {input.graph} {input.fasta} {output.fasta}"

rule dump_assembly_contigs_cull_and_no_derep:
    output:
        fasta="{stemA}.ggcat-{unitig_source}.{stemB}.contigs-d{min_depth}.fn",
    wildcard_constraints:
        unitig_source=noperiod_wc,
        min_depth=integer_wc,
    input:
        graph="{stemA}.ggcat-{unitig_source}.{stemB}.sz",
        fasta="{stemA}.ggcat-{unitig_source}.fn",
    params:
        min_depth=lambda w: int(w.min_depth) / 100
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip dump_contigs --verbose --min-depth {params.min_depth} --no-derep {input.graph} {input.fasta} {output.fasta}"


rule dump_assembly_segments:
    output:
        segments="{stema}.ggcat-{unitig_source}.{stemb}.segments.tsv",
    wildcard_constraints:
        unitig_source=noperiod_wc,
    input:
        graph="{stema}.ggcat-{unitig_source}.{stemb}.sz",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip dump_segments --verbose {input.graph} {output.segments}"


rule dump_assembly_depth:
    output:
        depth="{stema}.ggcat-{unitig_source}.{stemb}.sequence_depth.nc",
    wildcard_constraints:
        unitig_source=noperiod_wc,
    input:
        graph="{stema}.ggcat-{unitig_source}.{stemb}.sz",
        depth="{stema}.ggcat-{unitig_source}.unitig_depth.nc",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip dump_depth --verbose {input.graph} {input.depth} {output.depth}"


# rule dump_assembly_results:
#     output:
#         fasta="{stema}.ggcat-{unitig_source}.{stemb}.fn",
#         depth="{stema}.ggcat-{unitig_source}.{stemb}.sequence_depth.nc",
#         segments="{stema}.ggcat-{unitig_source}.{stemb}.segments.tsv",
#     wildcard_constraints:
#         unitig_source=noperiod_wc,
#     input:
#         graph="{stema}.ggcat-{unitig_source}.{stemb}.sz",
#         fasta="{stema}.ggcat-{unitig_source}.fn",
#         depth="{stema}.ggcat-{unitig_source}.unitig_depth.nc",
#     conda:
#         "conda/strainzip.yaml"
#     shell:
#         "strainzip dump --verbose {input.graph} {input.fasta} {input.depth} {output.segments} {output.depth} {output.fasta}"


rule megahit_assemble_single_k:
    output:
        dir=directory("data/group/{group}/r.{stem}.megahit-k{ksize}.d"),
        fasta="data/group/{group}/r.{stem}.megahit-k{ksize}.fn",
    wildcard_constraints:
        ksize=integer_wc,
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        r1=lambda w: ",".join(
            [
                f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
        r2=lambda w: ",".join(
            [
                f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
    resources:
        memory_flags="--memory 0.5",
    threads: 48
    shell:
        """
        megahit -t {threads} --k-list {wildcards.ksize} {params.memory_flags} -o {output.dir} -1 {params.r1} -2 {params.r2}
        ln {output.dir}/final.contigs.fa {output.fasta}
        """


rule megahit_assemble_k_series:
    output:
        dir=directory("data/group/{group}/r.{stem}.megahit-full-k{ksize}.d"),
        fasta="data/group/{group}/r.{stem}.megahit-full-k{ksize}.fn",
        fasta_unfiltered="data/group/{group}/r.{stem}.megahit-full-k{ksize}-unfiltered.fn",
    wildcard_constraints:
        ksize=integer_wc,
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        r1=lambda w: ",".join(
            [
                f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
        r2=lambda w: ",".join(
            [
                f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
        kmin=31,
        kmax=lambda w: w.ksize,
        kstep=20,
    threads: 48
    resources:
        memory_flags="--memory 0.5",
    shell:
        """
        megahit -t {threads} --k-min {params.kmin} --k-max {params.kmax} --k-step {params.kstep} \
                {resources.memory_flags} -o {output.dir} -1 {params.r1} -2 {params.r2}
        ln {output.dir}/final.contigs.fa {output.fasta}
        ln {output.dir}/intermediate_contigs/k{params.kmax}.contigs.fa {output.fasta_unfiltered}
        """


rule compile_ref_genome_mapping:
    output:
        "data/group/{group}/reference_genome_quast_contigs.tsv",
    input:
        script="scripts/compile_ref_genome_mapping_for_quast.py",
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    params:
        ref_args=lambda w: [
            f"{genome}=data/genome/{genome}.fn"
            for genome in config["genome_group"][w.group]
        ],
    shell:
        "{input.script} {params.ref_args} > {output}"

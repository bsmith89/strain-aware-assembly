config["genome"] = pd.read_table("meta/genome.tsv", dtype=str, index_col=["genome_id"])
config["genome_group"]["xjin"] = pd.read_table(
    "meta/genome.tsv", dtype=str, index_col=["genome_id"]
).index


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


rule link_project_reference_genome_no_species:
    output:
        "data/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


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


# NOTE: --hard-min 0 and --share-min 1 always
# or else some 0s will be censoring instead of actually
# missing. (That's because counts below --hard-min are discarded
# from individual samples and counts below --soft-min can be
# discarded from individual samples if they're solid in fewer than
# --share-min samples.
# With --hard-min 0 and --share-min 1, counts for a kmer are either
# kept entirely or discarded entirely.
rule run_kmtricks_pipeline:
    output:
        directory("{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d"),
    input:
        "{stem}.kmtricks_input.txt",
    wildcard_constraints:
        ksize=integer_wc,
        mincount=integer_wc,
        recurrence=integer_wc,
    params:
        ksize=lambda w: int(w.ksize),
        mincount=lambda w: int(w.mincount),
        recurrence=lambda w: int(w.recurrence),
        num_partitions=32,
    conda:
        "conda/kmtricks.yaml"
    threads: 24
    shell:
        """
        workdir={output}.tmp
        kmtricks pipeline --until superk --run-dir $workdir \
                --threads {threads} --verbose info \
                --kmer-size {params.ksize} \
                --nb-partitions {params.num_partitions} \
                --file {input}


        cut -d' ' -f1 {input} | xargs -n1 -I % \
            kmtricks count --run-dir $workdir \
                --threads {threads} --verbose info \
                --mode kmer \
                --hard-min 0 \
                --id %

        seq 0 31 | xargs -n1 -I % -P {threads} \
            kmtricks merge --run-dir $workdir \
                --threads 1 --verbose info \
                --mode kmer:count:text \
                --share-min 1 \
                --soft-min {params.mincount} \
                --recurrence-min {params.recurrence} \
                --partition-id %

        mv $workdir/matrices {output}
        mv $workdir {output}/workdir
        """


rule load_kmtricks_output_to_sqlite:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    input:
        script="scripts/sample_list_to_sqlite_table_script.py",
        counts="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
        sample_list="{stem}.kmtricks_input.txt",
    shell:
        """
        tmpdb={output}.tmp
        {input.script} {input.sample_list} {wildcards.ksize} | sqlite3 $tmpdb
        sort -m {input.counts}/matrix_*.count.txt | tqdm --unit-scale 1 | sqlite3 -separator ' ' $tmpdb '.import /dev/stdin count_'
        mv $tmpdb {output}
        """


rule run_ggcat_on_kmtricks_kmers:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat.fn",
    input:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
    params:
        quality_string=lambda w: "#",  # NOTE: Not correctly formatted?: * int(w.ksize)
    container:
        config["container"]["ggcat"]
    threads: 36
    shell:
        """
        fifo_dir=$(mktemp -d)
        for file in {input}/*.count.txt;
        do
            fq=$fifo_dir/$(basename $file).fifo.fq
            echo making fifo $fq
            mkfifo $fq
            awk '{{print $1}}' $file | sed 's:^\(.*\)$:@\\n\\1\\n+\\n{params.quality_string}:' > $fq &
            echo fifo $fq running
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $fifo_dir/*
        rm -r $fifo_dir
        """


rule calculate_mean_unitig_depths_across_samples:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat.unitig_depth.nc",
    input:
        fasta="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat.fn",
        db="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    conda:
        "conda/strainzip.yaml"
    threads: 12
    shell:
        """
        strainzip depth --preload -p {threads} {input.fasta} {input.db} {wildcards.ksize} {output}
        """


rule load_ggcat_fasta_to_sz:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat.sz",
    input:
        fasta="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat.fn",
        depth="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat.unitig_depth.nc",
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip load --verbose {wildcards.ksize} {input.fasta} {input.depth} {output}
        """


rule depth_smooth:
    output:
        "{stem}.smoothed.sz",
    input:
        "{stem}.sz",
    conda:
        "conda/strainzip.yaml"
    threads: 36
    shell:
        """
        strainzip smooth --verbose -p {threads} {input} {output}
        """


rule trim_tips:
    output:
        "{stem}.notips.sz",
    input:
        "{stem}.sz",
    conda:
        "conda/strainzip.yaml"
    threads: 36
    shell:
        "strainzip trim -p {threads} --verbose {input} {output}"


rule trim_tips_unpressed:
    output:
        "{stem}.notips-unpressed.sz",
    input:
        "{stem}.sz",
    conda:
        "conda/strainzip.yaml"
    threads: 36
    shell:
        "strainzip trim -p {threads} --verbose --no-press {input} {output}"


rule deconvolve_junctions:
    output: "{stem}.deconvolve-{thresh}.sz"
    input: "{stem}.sz"
    conda: "conda/strainzip.yaml"
    threads: 36
    shell: "strainzip assemble -p {threads} --verbose {input} {wildcards.thresh} {output}"


rule extract_assembly_results:
    output:
        fasta="{stemA}.ggcat.{stemB}.fn",
        depth="{stemA}.ggcat.{stemB}.depth.tsv",
        segments="{stemA}.ggcat.{stemB}.segments.tsv",
    input:
        graph="{stemA}.ggcat.{stemB}.sz",
        fasta="{stemA}.ggcat.fn",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip extract --verbose {input.graph} {input.fasta} {output.segments} {output.depth} {output.fasta}"


rule quality_asses_assembly_against_one_ref:
    output:
        directory("{stem}.gquast-{genome}.d"),
    input:
        tigs="{stem}.fn",
        ref="data/genome/{genome}.fn",
    threads: 24
    params:
        min_tig_length=1000,
    conda:
        "conda/quast.yaml"
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_tig_length} -R {input.ref} --output-dir {output} {input.tigs}
        """


rule quality_asses_assembly_against_all_refs:
    output:
        directory("{stem}.quast-{group}.d"),
    input:
        tigs="{stem}.fn",
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    threads: 24
    params:
        min_tig_length=1000,
        refs=lambda w: ",".join(
            [f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]]
        ),
    conda:
        "conda/quast.yaml"
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_tig_length} -r {params.refs} --output-dir {output} {input.tigs}
        """


# rule convert_bcalm_to_gfa:
#     output:
#         "{stem}.bcalm-k{ksize}.gfa",
#     input:
#         script="scripts/bcalm_to_gfa.py",
#         fn="{stem}.bcalm-k{ksize}.fn",
#     params:
#         ksize=lambda w: int(w.ksize),
#     shell:
#         "{input.script} {input.fn} {output} {params.ksize}"

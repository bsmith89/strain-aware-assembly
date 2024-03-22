rule start_shell_ggcat:
    container: config["container"]["ggcat"]
    shell:
        "bash"

rule start_conda_shell:
    output: "start_shell.{conda}"
    conda: "conda/{conda}.yaml"
    shell:
        "bash"

rule start_conda_ipython:
    output: "start_ipython.{conda}"
    conda: "conda/{conda}.yaml"
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
        dd("""
        cat <<EOF > {output}
        ecoli_mg1655 : data/genbank/ecoli.mg1655.fn
        ecoli_o121h19 : data/genbank/ecoli.o121h19.fn
        bdorei_dsm17855 : data/genbank/bdorei.dsm17855.fn
        EOF
        """)

rule construct_two_genome_input_table:
    output:
        "data/two_genomes.kmtricks_input.txt",
    input:
        ecoli_mg1655="data/genbank/ecoli.mg1655.fn",
        ecoli_o121h19="data/genbank/ecoli.o121h19.fn",
    shell:
        dd("""
        cat <<EOF > {output}
        ecoli_mg1655 : data/genbank/ecoli.mg1655.fn
        ecoli_o121h19 : data/genbank/ecoli.o121h19.fn
        EOF
        """)


rule run_kmtricks_pipeline:
    output:
        directory("{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d"),
    input:
        "{stem}.kmtricks_input.txt",
    params:
        ksize=lambda w: int(w.ksize),
        mincount=lambda w: int(w.mincount),
        recurrence=lambda w: int(w.recurrence),
    conda:
        "conda/kmtricks.yaml"
    threads: 24
    shell:
        """
        tmpdir={output}.tmp
        kmtricks pipeline \
                --kmer-size {params.ksize} \
                --hard-min 0 --share-min 1 --soft-min {params.mincount} --recurrence-min {params.recurrence} \
                --file {input} --run-dir $tmpdir \
                --mode kmer:count:text \
                --threads {threads}
        mv $tmpdir/matrices {output}
        rm -r $tmpdir
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
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    container:
        config["container"]["ggcat"]
    threads: 24
    shell:
        """
        fifo_dir=$(mktemp -d)
        for i in $(seq {threads});
        do
            fq=$fifo_dir/kmer.fifo-$i.fq
            echo making fifo $fq
            mkfifo $fq
            sqlite3 {input} "SELECT '@' || rowid || '|' || kmer || '|+|' || kmer AS rowid FROM count_ WHERE rowid % {threads} = ($i - 1);" \
                    | tr '|' '\\n' \
                > $fq &
            echo fifo $fq running
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $fifo_dir/kmer.fifo-*.fq
        rm -r $fifo_dir
        """


rule convert_bcalm_to_gfa:
    output:
        "{stem}.bcalm-k{ksize}.gfa",
    input:
        script="scripts/bcalm_to_gfa.py",
        fn="{stem}.bcalm-k{ksize}.fn",
    params:
        ksize=lambda w: int(w.ksize),
    shell:
        "{input.script} {input.fn} {output} {params.ksize}"

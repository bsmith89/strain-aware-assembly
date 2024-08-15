rule download_pgam_hmm_db:
    output:
        "raw/ref/hmm_PGAP.HMM.tgz",
    params:
        url="https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz",
    shell:
        curl_recipe


localrules:
    download_pgam_hmm_db,


rule extract_pgam_hmm_db:
    output:
        "ref/hmm/PGAP.hmm",
    input:
        "raw/ref/hmm_PGAP.HMM.tgz",
    shell:
        """
        tar -O -xzf {input} > {output}
        """


rule extract_single_tigrfam_model:
    output:
        "ref/hmm/TIGR{number}.hmm",
    input:
        "ref/hmm/PGAP.hmm",
    shell:
        """
        hmmfetch {input} TIGR{wildcards.number}.1 > {output}
        """
# NOTE (2024-06-19): I only accept TIGRXXXXX.1 here,
# even if there are later versions like .2


rule prokka_annotate_cluster:
    output:
        directory("data/{stem}.prokka.d"),
    input: "data/{stem}.fn"
    threads: 24
    conda: 'conda/prokka.yaml'
    shell:
        """
        prokka --force --cpus {threads} {input} \
                --outdir {output} \
                --prefix PROKKA \
                --locustag TODO \
                --cdsrnaolap
        """


rule run_resfinder:
    output:
        dir=directory("{stem}.resfinder.d")
    input:
        fasta="{stem}.fn",
        db="ref/resfinder_db",
    conda:
        "conda/resfinder.yaml"
    shell:
        "run_resfinder.py --acquired --db_path_res {input.db} --inputfasta {input.fasta} -o {output.dir}"


rule download_genomad_db:
    output:
        "ref/genomad_db"
    conda:
        "conda/genomad.yaml"
    shell:
        "genomad download-database $(dirname {output})"


rule run_genomad:
    output: directory("{stem}.genomad.d")  # FIXME
    input:
        db="ref/genomad_db",
        fasta="{stem}.fn",
    params:
        dir="{stem}.genomad.d"
    conda:
        "conda/genomad.yaml"
    threads: 48
    shell:
        """
        tmpdir=$(mktemp -d)
        ln -s $(realpath {input.fasta}) $tmpdir/contigs
        genomad end-to-end --threads {threads} --restart --verbose $tmpdir/contigs {params.dir} {input.db}
        rm -r $tmpdir
        """


rule download_mmseqs_uniref50_db:
    output:
        directory("ref/mmseqs_uniref50/"),
    conda:
        "conda/mmseqs.yaml"
    params:
        stem="ref/mmseqs_uniref50/db"
    threads: 4
    shell:
        "mkdir -p {output} && mmseqs databases --threads {threads} UniRef50 {params.stem} $TMPDIR"

# rule run_mmseqs_taxonomy:
#     output:


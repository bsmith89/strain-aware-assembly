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


# rule compare_crispr_annotations_minced_megahit_and_strainzip:
#     output:
#         sz=" data/group/watson_donor_a/r.proc.ggcat-k111-withmegahit2-droptips.notips-2.smoothed-6.unzip-norm-10-10.cctk-compare.d/MINCED_OUT/strainzip_minced_out.txt",
#         mh=" data/group/watson_donor_a/r.proc.ggcat-k111-withmegahit2-droptips.notips-2.smoothed-6.unzip-norm-10-10.cctk-compare.d/MINCED_OUT/megahit_minced_out.txt",
#     input:
#         sz=" data/group/watson_donor_a/r.proc.ggcat-k111-withmegahit2-droptips.notips-2.smoothed-6.unzip-norm-10-10.fn",
#         mh=" data/group/watson_donor_a/r.proc.megahit-full-k111.fn",
#     params:
#         cctkdir="data/group/watson_donor_a/r.proc.ggcat-k111-withmegahit2-droptips.notips-2.smoothed-6.unzip-norm-10-10.cctk-compare.d"
#     conda:
#         "conda/cctk.yaml"
#     shell:
#         """
#         tmpdir=$(mktemp -d) && echo $tmpdir
#         ln -s $(realpath {input.sz}) $tmpdir/strainzip.fn
#         ln -s $(realpath {input.mh}) $tmpdir/megahit.fn
#         cctk minced --run-minced -o {params.cctkdir} -i $tmpdir
#         rm -r $tmpdir
#         """


rule run_crispr_array_annotation_minced:
    output:
        minced="{stem}.cctk.d/MINCED_OUT/contigs_minced_out.txt",
    input:
        "{stem}.fn"
    params:
        cctkdir="{stem}.cctk.d"
    conda:
        "conda/cctk.yaml"
    shell:
        """
        tmpdir=$(mktemp -d) && echo $tmpdir
        ln -s $(realpath {input}) $tmpdir/contigs.fn
        cctk minced --run-minced -o {params.cctkdir} -i $tmpdir
        rm -r $tmpdir
        """


rule run_crispr_array_annotation_minced_postprocessing:
    output:
        processed=directory("{stem}.cctk.d/PROCESSED"),
    input:
        "{stem}.cctk.d/MINCED_OUT/contigs_minced_out.txt"
    params:
        cctkdir="{stem}.cctk.d",
        snp_thresh=0,
        min_shared=0,
    conda:
        "conda/cctk.yaml"
    shell:
        """
        cctk minced -o {params.cctkdir} -p --snp-thresh {params.snp_thresh} --min-shared {params.min_shared}
        """


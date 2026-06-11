rule download_pgam_hmm_db:
    output:
        "raw/ref/hmm_PGAP.HMM.tgz",
    params:
        url="https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz",
    shell:
        curl_recipe


localrules:
    download_pgam_hmm_db,


# rule extract_pgam_hmm_db:
#     output:
#         "ref/hmm/PGAP.hmm",
#     input:
#         "raw/ref/hmm_PGAP.HMM.tgz",
#     shell:
#         """
#         tar -O -xzf {input} > {output}
#         """
#
#
# rule extract_single_tigrfam_model:
#     output:
#         "ref/hmm/TIGR{number}.hmm",
#     input:
#         "ref/hmm/PGAP.hmm",
#     shell:
#         """
#         hmmfetch {input} TIGR{wildcards.number}.1 > {output}
#         """

rule extract_single_pgap_model:
    output:
        "ref/hmm/{model}_{ver}.hmm",
    input:
        "ref/hmm/PGAP.hmm",
    wildcard_constraints:
        model=single_param_wc,
        ver=single_param_wc,
    params:
        model_accession=lambda w: w.model + "." + w.ver
    shell:
        """
        hmmfetch {input} {params.model_accession} > {output}
        """


# NOTE (2024-06-19): I only accept TIGRXXXXX.1 here,
# even if there are later versions like .2


rule prokka_annotate_cluster:
    output:
        directory("data/{stem}.prokka.d"),
    input:
        "data/{stem}.fn",
    threads: 24
    conda:
        "conda/prokka.yaml"
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
        dir=directory("{stem}.resfinder.d"),
    input:
        fasta="{stem}.fn",
        db="ref/resfinder_db",
    conda:
        "conda/resfinder.yaml"
    shell:
        "run_resfinder.py --acquired --db_path_res {input.db} --inputfasta {input.fasta} -o {output.dir}"


rule download_genomad_db:
    output:
        "ref/genomad_db",
    conda:
        "conda/genomad.yaml"
    shell:
        "genomad download-database $(dirname {output})"


rule run_genomad:
    output:
        directory("{stem}.genomad.d"),  # FIXME
    input:
        db="ref/genomad_db",
        fasta="{stem}.fn",
    params:
        dir="{stem}.genomad.d",
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


rule download_mmseqs_db:
    output:
        "ref/mmseqs/{db}",
    threads: 12
    conda:
        "conda/mmseqs.yaml"
    shell:
        """
        mmseqs databases --threads {threads} {wildcards.db} {output} {resources.tmpdir}
        """

rule run_mmseqs_taxonomy:
    output:
        directory("{stem}.mmseqs-{db}-taxonomy.d"),
    input:
        db="ref/mmseqs/{db}",
        query="{stem}.fn",
    params:
        sensitivity=2.0,
    resources:
        disk_mb=10_000,
    conda:
        "conda/mmseqs.yaml"
    threads: 64
    shell: """
    mkdir -p {output}
    mmseqs easy-taxonomy --threads {threads} --disk-space-limit {resources.disk_mb}M -s {params.sensitivity} {input.query} {input.db} {output}/mmseqs {resources.tmpdir}
    """



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
        "{stem}.fn",
    params:
        cctkdir="{stem}.cctk.d",
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
        "{stem}.cctk.d/MINCED_OUT/contigs_minced_out.txt",
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

rule eggnog_mapper_translated_orfs:
    output:
        directory("{stem}.emapper.d"),
    input:
        fasta="{stem}.tran.fa",
        db="ref/eggnog_mapper_db",
    params:
        tax_scope="auto",
        sensmode="more-sensitive",
        mapper="diamond",
    conda:
        "conda/emapper.yaml"
    threads: 48
    resources:
        walltime_hr=240,
        mem_mb=20_000,
        pmem=20_000 // 48,
    shell:
        """
        tmpdir=$(mktemp -d)
        export EGGNOG_DATA_DIR={input.db}
        rm -rf {output}.temp
        mkdir -p {output}.temp
        emapper.py \
                -m {params.mapper} \
                -i {input.fasta} \
                --itype proteins \
                --sensmode {params.sensmode} \
                --go_evidence all \
                --dbmem \
                --tax_scope {params.tax_scope} \
                --temp_dir $tmpdir \
                --override \
                --cpu {threads} \
                --output_dir {output}.temp \
                --output 'proteins'
        rm -rf {output}
        mv {output}.temp {output}
        """

rule parse_strain_emapper_annotations_to_gene_x_unit:
    output:
        "{stem}.cds.emapper.gene_x_{unit}.tsv",
    input:
        script="scripts/parse_emapper_output_to_gene_x_{unit}.py",
        emapper="{stem}.cds.emapper.d/proteins.emapper.annotations",
    shell:
        "{input.script} {input.emapper} {output}"


rule gather_gene_seq_for_cog:
    output:
        nucl="{stem}.cds.{cog}.fn",
    wildcard_constraints:
        cog="COG[0-9]+"
    input:
        emapper_cog_annot="{stem}.cds.emapper.gene_x_cog.tsv",
        nucl="{stem}.cds.fn",
    conda:
        "conda/seqtk.yaml"
    shell:
        """
        seqtk subseq {input.nucl} <(grep '{wildcards.cog}' {input.emapper_cog_annot} | cut -f1) > {output.nucl}
        """

rule codonalign:
    output: "{stem}.codonalign.afn"
    input:
        prot="{stem}.tran.muscle.afa",
        nucl="{stem}.fn"
    conda:
        "conda/compbio_scripts.yaml"
    shell: "codonalign {input.prot} {input.nucl} > {output}"


rule infer_codon_phylogeny:
    output: "{stem}.codonalign.nucl.nwk",
    input: "{stem}.codonalign.afn"
    shell: "fasttree -nt < {input} > {output}"

rule tree_sort_afn:
    output: "{stem}.tree-sort.afn"
    input:
        script="scripts/get_ordered_leaves.py",
        tree="{stem}.nucl.nwk",
        seqs="{stem}.afn"
    conda:
        "conda/compbio_scripts.yaml"
    shell: "fetch_seqs --match-order <({input.script} {input.tree}) {input.seqs} > {output}"

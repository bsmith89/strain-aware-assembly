rule quality_asses_assembly_against_all_refs:
    output:
        dir=directory("{stem}.quast-{group}.d"),
        contig_to_genome="{stem}.quast-{group}.contig_to_genome.txt",
    input:
        tigs="{stem}.fn",
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    threads: 48
    params:
        min_tig_length=1000,
        min_identity=99,
        min_alignment=100,  # FIXME: Should really be ksize
        refs=lambda w: ",".join(
            [f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]]
        ),
    conda:
        "conda/quast.yaml"
    shell:
        """
        metaquast.py --silent --threads={threads} \
                --fragmented --min-contig {params.min_tig_length} \
                --min-alignment {params.min_alignment} --min-identity {params.min_identity} \
                -r {params.refs} \
                --output-dir {output.dir} \
                {input.tigs}
        cp {output.dir}/combined_reference/contigs_reports/alignments_*.tsv {output.contig_to_genome}
        # Reduce storage requirements:
        rm -rf {output.dir}/{{combined_reference/icarus_viewers,quast_corrected_input,runs_per_reference}}
        """


rule quality_asses_assembly_against_one_ref:
    output:
        dir=directory("{stem}.gquast-{genome}.d"),
    input:
        tigs="{stem}.fn",
        ref="data/genome/{genome}.fn",
    threads: 12
    params:
        min_tig_length=100,
        min_identity=95,
        min_alignment=50,  # FIXME: Should really be ksize
    conda:
        "conda/quast.yaml"
    shell:
        """
        quast --silent --threads={threads} \
                --fragmented --min-contig {params.min_tig_length} \
                --min-alignment {params.min_alignment} --min-identity {params.min_identity} \
                -r {input.ref} \
                --output-dir {output.dir} \
                {input.tigs}
        """


# rule quality_asses_assembly_against_one_ref:
#     output:
#         directory("{stem}.gquast-{genome}.d"),
#     input:
#         tigs="{stem}.fn",
#         ref="data/genome/{genome}.fn",
#     threads: 48
#     params:
#         min_tig_length=1000,
#     conda:
#         "conda/quast.yaml"
#     shell:
#         """
#         metaquast.py --threads={threads} --min-contig {params.min_tig_length} -R {input.ref} --output-dir {output} {input.tigs}
#         """
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


rule combine_reference_genomes:
    output:
        "data/group/{group}/all_refs.fn",
    input:
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    threads: 1
    shell:
        "cat {input} > {output}"


# NOTE (2024-11-20): This is one example of a benchmark_sequences files. I could also do
# other things like 16S or arbitrary fragments.
rule combine_reference_cds_as_benchmark_sequences:
    output:
        "data/group/{group}/all_ref_cds.benchmark_sequences.fn",
    input:
        refs=lambda w: [
            f"data/genome/{genome}.cds.fn"
            for genome in config["genome_group"][w.group]
        ],
    shell:
        "cat {input} > {output}"


rule combine_reference_tiles_as_benchmark_sequences:
    output:
        "data/group/{group}/all_ref_tiles-k{ksize}-o{overlap}.benchmark_sequences.fn",
    input:
        refs=lambda w: [
            f"data/genome/{genome}.norm.tiles-k{w.ksize}-o{w.overlap}.fn"
            for genome in config["genome_group"][w.group]
        ],
    shell:
        "cat {input} > {output}"


rule blastn_to_benchmark_sequences:
    output:
        "data/group/{group}/{stem}.bench-{benchmark_name}-blastn.tsv",
    input:
        query="data/group/{group}/{benchmark_name}.benchmark_sequences.fn",
        subject="data/group/{group}/{stem}.fn",
    shell:
        """
        blastn -subject {input.subject} -query {input.query} -outfmt 6 -max_target_seqs 1000000000 -perc_identity 50 > {output}
        """


rule collect_assembly_and_annotations:
    output:
        "data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.ASSEMBLY_DETAILS.flag",
    input:
        final="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.sz",
        quast="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.quast-{group}.d",
        cctk="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.cctk.d/PROCESSED",
        clust_meta="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.clust-e50-n20000-d20.meta.tsv",
        clust_vertex="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.clust-e50-n20000-d20.vertex.tsv",
        clust_segment="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.clust-e50-n20000-d20.segment.tsv",
        tigr02013="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.cds.tran.hmmer-TIGR02013-ga.tsv",
        resfinder="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.resfinder.d",
        genomad="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.genomad.d",
        unpressed="data/group/{group}/r.proc.{stemA}.notips-2-unpressed.sz",
    shell:
        "echo {input} > {output}"

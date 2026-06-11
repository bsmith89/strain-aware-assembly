rule quality_asses_assembly_against_group_refs:
    output:
        dir=directory("{stem}.quast-{group}.d"),
        contig_to_genome="{stem}.quast-{group}.contig_to_genome.txt",
    input:
        tigs="{stem}.fn",
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    threads:
        lambda w, input: len(input.refs)
    params:
        min_tig_length=1000,
        min_identity=95,
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

rule quality_statistics_assembly_against_group_refs_at_fixed_identity:
    output:
        dir=directory("{stem}.quast-{group}-i{ident}.d"),
    input:
        tigs="{stem}.fn",
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    threads: 4
    params:
        min_tig_length=200,
        min_identity=lambda w: int(w.ident) / 100,
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
                --strict-NA \
                -r {params.refs} \
                --space-efficient --no-plots --no-html --no-icarus --no-snps --fast \
                --output-dir {output.dir} \
                {input.tigs}
        """

ruleorder: quality_statistics_assembly_against_group_refs_at_fixed_identity > quality_asses_assembly_against_group_refs



rule quality_asses_assembly_against_one_ref:
    output:
        dir=directory("{stem}.gquast-{genome}.d"),
    input:
        tigs="{stem}.fn",
        ref="data/genome/{genome}.fn",
    threads: 2
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


rule blastn_benchmark_sequences_as_query:
    output:
        "data/group/{group}/{stem}.bench-{benchmark_name}-blastn.tsv",
    input:
        query="data/group/{group}/{benchmark_name}.benchmark_sequences.fn",
        subject="data/group/{group}/{stem}.fn",
    threads: 24
    shell:
        """
        cat {input.query} | parallel --gnu --plain -j {threads} --recstart '>' --pipe \
                blastn -subject {input.subject} -query - -outfmt 6 -max_target_seqs 1000000000 -perc_identity 50 \
        > {output}
        """

use rule blastn_benchmark_sequences_as_query as blastn_refs_as_subject with:
    output:
        "data/group/{group}/{stem}.all_refs-blastn.tsv",
    input:
        query="data/group/{group}/{stem}.fn",
        subject="data/group/{group}/all_refs.fn",


rule collect_assembly_and_annotations:
    output:
        "data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.ASSEMBLY_DETAILS.flag",
    input:
        final="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.sz",
        # quast0="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.quast-{group}.d",
        strainzip_quast_range=[f"data/group/{{group}}/r.proc.{{stemA}}.notips-2.{{stemB}}.contigs-d0.quast-{{group}}-i{ident}.d" for ident in ['10000', '9990', '9950', '9900', '9800']],
        # resfinder="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.resfinder.d",
        # genomad="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.genomad.d",
        # cctk="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.cctk.d/PROCESSED",
        # clust_meta="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.clust-e50-n20000-d20.meta.tsv",
        # clust_vertex="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.clust-e50-n20000-d20.vertex.tsv",
        # clust_segment="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.clust-e50-n20000-d20.segment.tsv",
        strainzip_tigr02013="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.cds.tran.hmmer-TIGR02013_1-ga.tsv",
        strainzip_tigr00952="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.cds.tran.hmmer-TIGR00952_1-ga.tsv",
        strainzip_tigr01063="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.cds.tran.hmmer-TIGR01063_1-ga.tsv",
        unpressed="data/group/{group}/r.proc.{stemA}.notips-2-unpressed.sz",
        megahit_quast_range=[f"data/group/{{group}}/r.proc.megahit-full-k111.quast-{{group}}-i{ident}.d" for ident in ['10000', '9990', '9950', '9900', '9800']],
        megahit_gene_depth="data/group/{group}/r.proc.megahit-full-k111.prodigal.gene_depth.tsv",
        megahit_gene_bed="data/group/{group}/r.proc.megahit-full-k111.prodigal.bed",
        megahit_contig_depth="data/group/{group}/r.proc.megahit-full-k111.contig_depth.tsv",
        megahit_tigr02013="data/group/{group}/r.proc.megahit-full-k111.cds.tran.hmmer-TIGR02013_1-ga.tsv",
        megahit_tigr00952="data/group/{group}/r.proc.megahit-full-k111.cds.tran.hmmer-TIGR00952_1-ga.tsv",
        megahit_tigr01063="data/group/{group}/r.proc.megahit-full-k111.cds.tran.hmmer-TIGR01063_1-ga.tsv",
        strainzip_gene_bed="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.prodigal.bed",
        tile10k_bench_nlength="data/group/{group}/all_ref_tiles-k10000-o5000.benchmark_sequences.nlength.tsv",
        tile10k_bench_matching="data/group/{group}/all_refs.bench-all_ref_tiles-k10000-o5000-blastn.tsv",
        tile10k_bench_strainzip="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.bench-all_ref_tiles-k10000-o5000-blastn.tsv",
        tile10k_bench_megahit="data/group/{group}/r.proc.megahit-full-k111.bench-all_ref_tiles-k10000-o5000-blastn.tsv",
        strainzip_to_refs="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.all_refs-blastn.tsv",
        strainzip_nlength="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.nlength.tsv",
        megahit_to_refs="data/group/{group}/r.proc.megahit-full-k111.all_refs-blastn.tsv",
        megahit_nlength="data/group/{group}/r.proc.megahit-full-k111.nlength.tsv",
        gene_bench_nlength="data/group/{group}/all_ref_cds.benchmark_sequences.nlength.tsv",
        gene_bench_matching="data/group/{group}/all_refs.bench-all_ref_cds-blastn.tsv",
        gene_bench_strainzip="data/group/{group}/r.proc.{stemA}.notips-2.{stemB}.contigs-d0.bench-all_ref_cds-blastn.tsv",
        gene_bench_megahit="data/group/{group}/r.proc.megahit-full-k111.bench-all_ref_cds-blastn.tsv",
    shell:
        "echo {input} > {output}"

rule bowtie_index_build_megahit_coassembly:
    output:
        bt2_1="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.1.bt2",
        bt2_2="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.2.bt2",
        bt2_3="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.3.bt2",
        bt2_4="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.4.bt2",
        bt2_r1="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.rev.1.bt2",
        bt2_r2="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.rev.2.bt2",
    input:
        "data/group/{group}/r.proc.megahit-full-k111.fn",
    params:
        bowtie2_index="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs",
    threads: 24
    resources:
        mem_mb=int(100e3),
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.bowtie2_index}
        """


rule bowtie2_map_reads_to_megahit_assembly:
    output:
        "data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.bam",
    input:
        r1="data/reads/{mgen}/r1.proc.fq.gz",
        r2="data/reads/{mgen}/r2.proc.fq.gz",
        bt2_1="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.1.bt2",
        bt2_2="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.2.bt2",
        bt2_3="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.3.bt2",
        bt2_4="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.4.bt2",
        bt2_r1="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.rev.1.bt2",
        bt2_r2="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs.rev.2.bt2",
    params:
        bowtie2_index="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs",
    conda:
        "conda/bowtie2.yaml"
    threads: 4
    shell:
        """
        tmp={output}.tmp
        echo "Writing temporary bamfile to $tmp for {output}"
        bowtie2 --threads {threads} \
                -x {params.bowtie2_index} \
                -1 {input.r1} -2 {input.r2} \
            | samtools view -b - \
            | samtools sort -o $tmp
        mv $tmp {output}
        """

rule tally_position_depth_from_bam:
    output:
        pos_depth="data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.position_depth.tsv",
    input:
        bam="data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.bam",
    shell:
        "samtools depth -a -s {input.bam} > {output}"


rule calculate_gene_mean_mapping_depth:
    output:
        "data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.prodigal.gene_depth.tsv",
    input:
        script="scripts/calculate_gene_mean_mapping_depth.sh",
        gff="data/group/{group}/r.proc.megahit-full-k111.prodigal.gff",
        pos_depth="data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.position_depth.tsv",
    threads: 1
    shell:
        """
        {input.script} {input.gff} {input.pos_depth} > {output}
        """

rule mean_depth_by_contig:
    output:
        "data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.contig_depth.tsv",
    input:
        pos_depth="data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.position_depth.tsv",
    threads: 1
    shell:
        """
            cat {input.pos_depth} | awk -v OFS='\t' '\\
                BEGIN {{gene_id="__START__"; depth_tally=0; pos_tally=0}} \\
                $1==gene_id {{depth_tally+=$3; pos_tally+=1}} \\
                $1!=gene_id {{print gene_id, 1.0 * depth_tally / pos_tally; gene_id=$1; depth_tally=0; pos_tally=0}} \\
                END {{print gene_id, 1.0 * depth_tally / pos_tally}} \\
                ' \
            | sed '1,1d' > {output}
        """

rule combine_gene_depth_across_group_samples:
    output: "data/group/{group}/r.proc.megahit-full-k111.prodigal.gene_depth.tsv"
    input:
        lambda w: [
            f"data/group/{w.group}/reads/{mgen}/r.proc.megahit-full-k111.prodigal.gene_depth.tsv"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        pattern="data/group/{group}/reads/$sample/r.proc.megahit-full-k111.prodigal.gene_depth.tsv"
    shell:
        """
        for sample in {params.mgen_list}
        do
            sample_depth={params.pattern}
            echo $sample >&2
            awk -v OFS='\t' -v sample=$sample '{{print sample,$1,$2}}' $sample_depth
        done > {output}
        """

rule combine_contig_depth_across_group_samples:
    output: "data/group/{group}/r.proc.megahit-full-k111.contig_depth.tsv"
    input:
        lambda w: [
            f"data/group/{w.group}/reads/{mgen}/r.proc.megahit-full-k111.contig_depth.tsv"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        pattern="data/group/{group}/reads/$sample/r.proc.megahit-full-k111.contig_depth.tsv"
    shell:
        """
        for sample in {params.mgen_list}
        do
            sample_depth={params.pattern}
            echo $sample >&2
            awk -v OFS='\t' -v sample=$sample '{{print sample,$1,$2}}' $sample_depth
        done > {output}
        """

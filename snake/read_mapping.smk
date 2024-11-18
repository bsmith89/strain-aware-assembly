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
        bowtie2_index="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs"
    threads: 24
    resources:
        mem_mb=int(100e3),
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.bowtie2_index}
        """

rule bowtie2_map_reads_to_megahit_assembly:
    output: 'data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.bam'
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
        bowtie2_index="data/group/{group}/r.proc.megahit-full-k111.bowtie2_index.d/contigs"
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


rule count_mapping_depth_along_contigs:
    output: 'data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.position_depth.tsv'
    input:
        bam="data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.bam",
    conda: "conda/bowtie2.yaml"
    threads: 1
    shell:
        """
        samtools depth --threads {threads} {input.bam} > {output}
        """


rule calculate_gene_mean_mapping_depth:
    output: 'data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.prodigal.gene_depth.tsv'
    input:
        script="scripts/calculate_gene_mean_mapping_depth.py",
        gff="data/group/{group}/r.proc.megahit-full-k111.prodigal.gff",
        pos_depth="data/group/{group}/reads/{mgen}/r.proc.megahit-full-k111.position_depth.tsv",
    threads: 1
    shell:
        """
        {input.script} {input.gff} {input.pos_depth} {output}
        """

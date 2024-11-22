rule download_sra_raw_fastq:
    output:
        r1="raw/sra/{srr}_1.fastq",
        r2="raw/sra/{srr}_2.fastq",
    conda:
        "conda/sra_tools.yaml"
    shell:
        "fasterq-dump --outdir raw/sra {wildcards.srr}"


rule gzip_sra_fastq:
    output:
        "raw/sra/{stem}.fq.gz",
    input:
        "raw/sra/{stem}.fastq",
    shell:
        "gzip -c {input} > {output}"

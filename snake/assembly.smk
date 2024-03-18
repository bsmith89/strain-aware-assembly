use rule start_shell as start_shell_bcalm with:
    conda:
        "conda/bcalm.yaml"


use rule start_shell as start_shell_jfish with:
    conda:
        "conda/jellyfish.yaml"


# FIXME: Cannot get python bindings to work.
use rule install_jupyter_kernel_default as install_jupyter_kernel_jfish with:
    params:
        name='jfish'
    conda:
        "conda/jellyfish.yaml"

def genbank_genomic_ftp_url(accession, assembly):
    prefix=accession[:3]
    n1to3=accession[4:7]
    n4to6=accession[7:10]
    n7to9=accession[10:13]
    return f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{n1to3}/{n4to6}/{n7to9}/{accession}_{assembly}/{accession}_{assembly}_genomic.fna.gz'
# e.g. https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/512/915/GCA_000512915.1_ASM51291v1/GCA_000512915.1_ASM51291v1_genomic.fna.gz

rule download_genbank_genome:
    output: 'raw/genbank/{accession}_{assembly}.fn'
    params:
        url=lambda w: genbank_genomic_ftp_url(w.accession, w.assembly)
    shell: "curl '{params.url}' | zcat > {output}"
localrules: download_genbank_genome


def alias_genbank_genome_input(w):
    accession, assembly = config['genomes'].loc[(w.species, w.strain)][['genbank', 'assembly']]
    return f'raw/genbank/{accession}_{assembly}.fn'


rule alias_genbank_genome:
    output: 'data/genbank/{species}.{strain}.fn'
    input: alias_genbank_genome_input
    shell: alias_recipe


rule merge_ecoli_strains:
    output: "data/both_strains.fn"
    input:
        "data/genbank/ecoli.mg1655.fn",
        "data/genbank/ecoli.o121h19.fn",
    shell:
        "cat {input} > {output}"


rule merge_ecoli_and_bdorei:
    output: "data/both_species.fn"
    input:
        "data/genbank/ecoli.mg1655.fn",
        "data/genbank/bdorei.dsm17855.fn",
    shell:
        "cat {input} > {output}"


rule merge_two_ecoli_and_bdorei:
    output: "data/three_genomes.fn"
    input:
        "data/genbank/ecoli.mg1655.fn",
        "data/genbank/ecoli.o121h19.fn",
        "data/genbank/bdorei.dsm17855.fn",
    shell:
        "cat {input} > {output}"


rule run_bcalm:
    output: '{stem}.bcalm-k{ksize}.fn'
    input: '{stem}.fn'
    params:
        outprefix='{stem}.bcalm-k{ksize}',
        ksize=lambda w: int(w.ksize),
    conda: 'conda/bcalm.yaml'
    threads: 24
    shell:
        """
        bcalm \
            -nb-cores {threads} \
            -in {input} \
            -kmer-size {params.ksize} \
            -abundance-min 1 \
            -out {params.outprefix}
        mv {params.outprefix}.unitigs.fa {output}
        """


rule convert_bcalm_to_gfa:
    output: "{stem}.bcalm-k{ksize}.gfa"
    input:
        script='scripts/bcalm_to_gfa.py',
        fn="{stem}.bcalm-k{ksize}.fn",
    params:
        ksize=lambda w: int(w.ksize),
    shell: "{input.script} {input.fn} {output} {params.ksize}"


rule run_jellyfish_count:
    output: "{stem}.jfish-k{ksize}.jf"
    input: "{stem}.fn"
    params:
        ksize=lambda w: int(w.ksize),
    threads: 4
    conda: 'conda/jellyfish.yaml'
    shell:
        """
        jellyfish count \
            --size 100M \
            --mer-len={params.ksize} \
            --threads={threads} \
            --canonical \
            --lower-count=0 \
            --output={output} \
            {input}
        """

rule dump_jellyfish_kmer_counts:
    output: "{stem}.kmer-k{ksize}.counts.tsv"
    input: "{stem}.jfish-k{ksize}.jf"
    conda: "conda/jellyfish.yaml"
    shell:
        """
        jellyfish dump --column --tab --output={output} {input}
        """

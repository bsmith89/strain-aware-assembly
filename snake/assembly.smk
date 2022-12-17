use rule start_shell as start_shell_bcalm with:
    conda:
        "conda/bcalm.yaml"

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
    species_id = w.species_id
    strain = w.strain
    accession, assembly = config['genomes'].loc[(species_id, strain)][['genbank', 'assembly']]
    return f'raw/genbank/{accession}_{assembly}.fn'


rule alias_genbank_genome:
    output: 'data/genbank/{species_id}.{strain}.fn'
    input: alias_genbank_genome_input
    shell: alias_recipe


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

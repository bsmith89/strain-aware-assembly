use rule start_shell as start_shell_bcalm with:
    conda:
        "conda/bcalm.yaml"


# TODO: Figure out how to link the dna_jellyfish package
# into my conda environment.
# NOTE: While I cannot get the conda environment building to automatically
# add the (swig) python bindings for jellyfish, I was able to figure out how
# to add them manually myself for a new conda environment being run inside of
# the jellyfish image:
# docker://bsmith89/jellyfish:435ac5c9499c6e4c92936884c5e549cd73816fdc
# To do this, you need to link three different files originally
# installed by pip into the base python for the Docker container.
# These files are all found in /opt/conda/lib/python3.9/site-packages/ :
#   dna_jellyfish-0.0.1.dist-info
#   dna_jellyfish.py
#   _dna_jellyfish.cpython-39-x86_64-linux-gnu.so
# All three should be linked into $CONDA_PREFIX/lib/python3.9/site-packages
# Unfortunately, I was not able to get this to work with *.post-deploy.sh
# because $CONDA_PREFIX was set to my outermost snakemake environment.
# (I think this is a bug.)
use rule start_shell as start_shell_jfish with:
    conda:
        "conda/jfish.yaml"
    container:
        config["container"]["jfish"]


# NOTE: This is mostly intended for testing jellyfish-dependent
# work while I still haven't figured out how to install it
# normally to a conda environment.
use rule start_shell as start_shell_jfish_base with:
    container:
        config["container"]["jfish"]


# NOTE See start_shell_jfish for how to get python bindings to work
# after building the conda environment.
# FIXME: Still cannot get bindings to work, even if they work fine
# from inside the singularity container.
use rule install_jupyter_kernel_native as install_jupyter_kernel_jfish with:
    params:
        name="jfish",
    conda:
        "conda/jfish.yaml"
    container:
        config["container"]["jfish"]


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
    output:
        "{stem}.kmer-k{ksize}.counts.tsv",
    input:
        "{stem}.jfish-k{ksize}.jf",
    container:
        config["container"]["jfish"]
    shell:
        """
        jellyfish dump --column --tab --output={output} {input}
        """

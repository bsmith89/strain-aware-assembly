config["genome"] = pd.read_table("meta/genome.tsv", dtype=str, index_col=["genome_id"])
config["genome_group"]["xjin"] = pd.read_table(
    "meta/genome.tsv", dtype=str, index_col=["genome_id"]
).index


rule start_shell_ggcat:
    container:
        config["container"]["ggcat"]
    shell:
        "bash"


rule start_conda_shell:
    output:
        "start_shell.{conda}",
    conda:
        "conda/{conda}.yaml"
    shell:
        "bash"


rule start_conda_ipython:
    output:
        "start_ipython.{conda}",
    conda:
        "conda/{conda}.yaml"
    shell:
        "ipython"


#
#
# use rule start_shell as start_shell_ggcat with:
#     container:
#         config["container"]["ggcat"]
#
#
# use rule start_shell as start_shell_kmtricks with:
#     conda:
#         "conda/kmtricks.yaml"
#


def genbank_genomic_ftp_url(accession, assembly):
    prefix = accession[:3]
    n1to3 = accession[4:7]
    n4to6 = accession[7:10]
    n7to9 = accession[10:13]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{n1to3}/{n4to6}/{n7to9}/{accession}_{assembly}/{accession}_{assembly}_genomic.fna.gz"


# e.g. https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/512/915/GCA_000512915.1_ASM51291v1/GCA_000512915.1_ASM51291v1_genomic.fna.gz


rule download_genbank_genome:
    output:
        "raw/genbank/{accession}_{assembly}.fn",
    params:
        url=lambda w: genbank_genomic_ftp_url(w.accession, w.assembly),
    shell:
        "curl '{params.url}' | zcat > {output}"


localrules:
    download_genbank_genome,


def alias_genbank_genome_input(w):
    accession, assembly = config["genomes"].loc[(w.species, w.strain)][
        ["genbank", "assembly"]
    ]
    return f"raw/genbank/{accession}_{assembly}.fn"


rule alias_genbank_genome:
    output:
        "data/genbank/{species}.{strain}.fn",
    input:
        alias_genbank_genome_input,
    shell:
        alias_recipe


rule link_project_reference_genome_no_species:
    output:
        "data/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


rule construct_mgen_group_input_table:
    output:
        "data/group/{group}/r.proc.kmtricks_input.txt",
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.proc.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.proc.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        entry_pattern="$mgen : data/reads/$mgen/r1.proc.fq.gz ; data/reads/$mgen/r2.proc.fq.gz",
        mgen_list=lambda w: config["mgen_group"][w.group],
    shell:
        """
        for mgen in {params.mgen_list}; do echo "{params.entry_pattern}"; done > {output}
        """


rule construct_three_genome_input_table:
    output:
        "data/three_genomes.kmtricks_input.txt",
    input:
        ecoli_mg1655="data/genbank/ecoli.mg1655.fn",
        ecoli_o121h19="data/genbank/ecoli.o121h19.fn",
        bdorei_dsm17855="data/genbank/bdorei.dsm17855.fn",
    shell:
        dd(
            """
        cat <<EOF > {output}
        ecoli_mg1655 : data/genbank/ecoli.mg1655.fn
        ecoli_o121h19 : data/genbank/ecoli.o121h19.fn
        bdorei_dsm17855 : data/genbank/bdorei.dsm17855.fn
        EOF
        """
        )


rule construct_two_genome_input_table:
    output:
        "data/two_genomes.kmtricks_input.txt",
    input:
        ecoli_mg1655="data/genbank/ecoli.mg1655.fn",
        ecoli_o121h19="data/genbank/ecoli.o121h19.fn",
    shell:
        dd(
            """
        cat <<EOF > {output}
        ecoli_mg1655 : data/genbank/ecoli.mg1655.fn
        ecoli_o121h19 : data/genbank/ecoli.o121h19.fn
        EOF
        """
        )


rule run_kmtricks_repart:
    output:
        repart_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/repart.flag",
    input:
        "{stem}.kmtricks_input.txt",
    params:
        workdir="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
        ksize=lambda w: int(w.ksize),
        num_partitions=32,  # This number must match the number used in compiling all the merge partitions.
    conda:
        "conda/kmtricks.yaml"
    threads: 8
    shell:
        """
        rm -r {params.workdir}  # Snakemake auto-builds this directory, so I need to delete it for kmtricks to run.
        kmtricks repart --run-dir {params.workdir} \
                --threads {threads} --verbose info \
                --kmer-size {params.ksize} \
                --nb-partitions {params.num_partitions} \
                --file {input}
        mkdir {params.workdir}/checkpoints
        touch {output.repart_flag}
        """


rule run_kmtricks_superk:
    output:
        superk_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/{mgen}.superk.flag",
    input:
        repart_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/repart.flag",
    params:
        workdir="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
    conda:
        "conda/kmtricks.yaml"
    threads: 1
    shell:
        """
        kmtricks superk --threads {threads} --verbose info --run-dir {params.workdir} --id {wildcards.mgen}
        touch {output.superk_flag}
        """


rule run_kmtricks_count:
    output:
        count_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/{mgen}.count.flag",
    input:
        superk_flag="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/checkpoints/{mgen}.superk.flag",
    params:
        workdir="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
    conda:
        "conda/kmtricks.yaml"
    threads: 8
    shell:
        """
        kmtricks count --run-dir {params.workdir} \
                --threads {threads} --verbose info \
                --mode kmer \
                --hard-min 0 \
                --id {wildcards.mgen}
        touch {output.count_flag}
        """


# NOTE: --hard-min 0 and --share-min 1 always
# or else some 0s will be censoring instead of actually
# missing. (That's because counts below --hard-min are discarded
# from individual samples and counts below --soft-min can be
# discarded from individual samples if they're solid in fewer than
# --share-min samples.
# With --hard-min 0 and --share-min 1, counts for a kmer are either
# kept entirely or discarded entirely.
rule run_kmtricks_merge:
    output:
        merge_flag="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d/matrices/matrix_{part}.count.txt",
    input:
        count_flags=lambda w: [
            f"data/group/{w.group}/{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/checkpoints/{mgen}.count.flag"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        workdir="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.d",
        mincount=lambda w: int(w.mincount),
        recurrence=lambda w: int(w.recurrence),
    conda:
        "conda/kmtricks.yaml"
    threads: 36
    shell:
        """
        kmtricks merge --run-dir {params.workdir} \
            --threads 1 --verbose info \
            --mode kmer:count:text \
            --share-min 1 \
            --soft-min {params.mincount} \
            --recurrence-min {params.recurrence} \
            --partition-id {wildcards.part}
        touch {output.merge_flag}
        """


rule load_kmtricks_output_to_sqlite:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    input:
        script="scripts/sample_list_to_sqlite_table_script.py",
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
        sample_list="{stem}.kmtricks_input.txt",
    shell:
        """
        tmpdb=$(mktemp) && echo $tmpdb
        {input.script} {input.sample_list} {wildcards.ksize} | sqlite3 $tmpdb
        sort -m {input.counts} | tqdm --unit-scale 1 | sqlite3 -separator ' ' $tmpdb '.import /dev/stdin count_'
        mv $tmpdb {output}
        """


rule run_ggcat_on_kmtricks_kmers:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-denovo.fn",
    input:
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
    container:
        config["container"]["ggcat"]
    threads: 36
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        for file in {input.counts}
        do
            fasta=$input_dir/$(basename $file).fifo.fa
            echo making fifo $fasta
            mkfifo $fasta
            awk -v OFS="\\n" '{{print ">"NR,$1}}' $file > $fasta &
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_ggcat_on_kmtricks_kmers_and_include_xjin_refs:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-withref.fn",
    input:
        counts=lambda w: [
            f"{w.stem}.kmtricks-k{w.ksize}-m{w.mincount}-r{w.recurrence}.d/matrices/matrix_{part}.count.txt"
            for part in range(32)
        ],
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"]["xjin"]
        ],
    container:
        config["container"]["ggcat"]
    threads: 36
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        for file in {input.refs}
        do
            ln -s $(realpath $file) $input_dir/$(basename $file).fa
        done

        for file in {input.counts}
        do
            fasta=$input_dir/$(basename $file).fifo.fa
            echo making fifo $fasta
            mkfifo $fasta
            awk -v OFS="\\n" '{{print ">"NR,$1}}' $file > $fasta &
        done
        sleep 2
        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule run_ggcat_on_xjin_refs:
    output:
        "data/xjin_ref.ggcat-k{ksize}-refonly.fn",
    input:
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"]["xjin"]
        ],
    container:
        config["container"]["ggcat"]
    threads: 36
    shell:
        """
        input_dir=$(mktemp -d)
        echo $input_dir

        for file in {input.refs}
        do
            ln -s $(realpath $file) $input_dir/$(basename $file).fa
        done

        ggcat build -s 1 -j {threads} -o {output} -k {wildcards.ksize} -e $input_dir/*
        rm -r $input_dir
        """


rule calculate_mean_unitig_depths_across_samples:
    output:
        "data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.unitig_depth.nc",
    wildcard_constraints:
        unitig_source=noperiod_wc,
    input:
        fasta="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.fn",
        db="data/group/{group}/{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.db",
    conda:
        "conda/strainzip.yaml"
    params:
        sample_list=lambda w: ",".join(config["mgen_group"][w.group]),
    threads: 36
    shell:
        """
        tmpdb=$(mktemp) && echo $tmpdb
        strainzip depth --verbose --preload --tmpdb $tmpdb -p {threads} {input.fasta} {input.db} {wildcards.ksize} {params.sample_list} {output}
        """


rule load_ggcat_with_depths_to_sz:
    output:
        "{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.sz",
    wildcard_constraints:
        unitig_source=single_param_wc,
    input:
        fasta="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.fn",
        depth="{stem}.kmtricks-k{ksize}-m{mincount}-r{recurrence}.ggcat-{unitig_source}.unitig_depth.nc",
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip load --verbose {wildcards.ksize} {input.fasta} {input.depth} {output}
        """


rule trim_tips:
    output:
        "{stem}.notips-{length}.sz",
    wildcard_constraints:
        length=single_param_wc,
    input:
        "{stem}.sz",
    params:
        length=lambda w: float(w.length),
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip trim --verbose -p {threads} --num-kmer-lengths {params.length} {input} {output}"


rule trim_tips_unpressed:
    output:
        "{stem}.notips-{length}-unpressed.sz",
    input:
        "{stem}.sz",
    params:
        length=lambda w: float(w.length),
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip trim --verbose -p {threads} --no-press --num-kmer-lengths {params.length} {input} {output}"


rule smooth_depths:
    output:
        "{stem}.smoothed-{eps}.sz",
    wildcard_constraints:
        eps=integer_wc,
    input:
        "{stem}.sz",
    params:
        eps=lambda w: 10**(-int(w.eps))
    conda:
        "conda/strainzip.yaml"
    threads: 36
    shell:
        """
        strainzip smooth --verbose -p {threads} --eps {params.eps} {input} {output}
        """


rule deconvolve_junctions:
    output:
        final="{stem}.deconvolve-{model}-{thresh}-{rounds}.sz",
    wildcard_constraints:
        model=single_param_wc,
        thresh=single_param_wc,
        rounds=single_param_wc,
    input:
        "{stem}.sz",
    log:
        checkpoint_dir=directory(
            "{stem}.deconvolve-{model}-{thresh}-{rounds}.checkpoints.d"
        ),
    params:
        model=lambda w: {
            "lognorm2": "OffsetLogNormal",
            "norm": "Normal",
            "lapl": "Laplace",
            "t5": "StudentsT --model-hyperparameters df=5",
            "huber": "Huber --model-hyperparameters delta=1",
        }[w.model],
        min_depth=1.0,
        score_thresh=lambda w: float(w.thresh),
        relative_error_thresh=0.1,
        absolute_error_thresh=1.0,
        max_rounds=lambda w: int(w.rounds),
        excess_thresh=0,
        completeness_thresh=1,
    conda:
        "conda/strainzip.yaml"
    threads: 36
    shell:
        """
        # FIXME: Figure out why setting environmental variables here is necessary.
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1 --xla_force_host_platform_device_count=8"
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREAD=1
        export NUM_INTER_THREADS=1
        export NUM_INTRA_THREADS=1

        mkdir -p {log.checkpoint_dir}

        strainzip deconvolve --verbose -p {threads} \
                --min-depth {params.min_depth} \
                --skip-extra-large --max-rounds {params.max_rounds} --model {params.model} \
                --score aic --score-thresh {params.score_thresh} \
                --relative-error-thresh {params.relative_error_thresh} \
                --absolute-error-thresh {params.absolute_error_thresh} \
                --excess-thresh {params.excess_thresh} \
                --completeness-thresh {params.completeness_thresh} \
                --checkpoint-dir {log.checkpoint_dir} \
                {input} {output.final}
        """


rule cluster_vertices:
    output:
        vertex="{stem}.clust-{thresh}.vertex.tsv",
        segment="{stem}.clust-{thresh}.segment.tsv",
        depth="{stem}.clust-{thresh}.depth.tsv",
        shared="{stem}.clust-{thresh}.shared.tsv",
        meta="{stem}.clust-{thresh}.meta.tsv",
    input:
        '{stem}.sz'
    params:
        num_preclust=10000,
        thresh=lambda w: int(w.thresh) / 1000,
        exponent=1/2,
    threads: 12
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip cluster --verbose \
                --num-preclust {params.num_preclust} --exponent {params.exponent} \
                {input} {params.thresh} {output.vertex} {output.segment} {output.depth} {output.shared} {output.meta}
        """

rule extract_unassembled_cluster_subgraph:
    output:
        graph="{stemA}.ggcat-{unitig_source}.{stemB}.clust-{thresh}.unassembled-c{clust}-r{radius}.sz",
    wildcard_constraints:
        radius=integer_wc,
    input:
        graph="{stemA}.ggcat-{unitig_source}.sz",
        segment="{stemA}.ggcat-{unitig_source}.{stemB}.clust-{thresh}.segment.tsv",
    params:
        clust=lambda w: int(w.clust),
        radius=lambda w: int(w.radius),
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip focus --verbose \
                --radius {params.radius} \
                --segments <(awk '$2 == {params.clust} {{print $1}}' {input.segment}) \
                {input.graph} \
                {output.graph}
        """

rule extract_assembled_cluster_subgraph:
    output:
        graph="{stem}.clust-{thresh}.assembled-c{clust}-r{radius}.sz",
    wildcard_constraints:
        radius=integer_wc,
    input:
        graph="{stem}.sz",
        vertex="{stem}.clust-{thresh}.vertex.tsv",
    params:
        clust=lambda w: int(w.clust),
        radius=lambda w: int(w.radius),
    conda:
        "conda/strainzip.yaml"
    shell:
        """
        strainzip focus --verbose \
                --radius {params.radius} \
                --vertices <(awk '$2 == {params.clust} {{print $1}}' {input.vertex}) \
                {input.graph} \
                {output.graph}
        """




rule extract_assembly_results:
    output:
        fasta="{stemA}.ggcat-{unitig_source}.{stemB}.fn",
        depth="{stemA}.ggcat-{unitig_source}.{stemB}.sequence_depth.nc",
        segments="{stemA}.ggcat-{unitig_source}.{stemB}.segments.tsv",
    wildcard_constraints:
        unitig_source=noperiod_wc,
    input:
        graph="{stemA}.ggcat-{unitig_source}.{stemB}.sz",
        fasta="{stemA}.ggcat-{unitig_source}.fn",
        depth="{stemA}.ggcat-{unitig_source}.unitig_depth.nc",
    conda:
        "conda/strainzip.yaml"
    shell:
        "strainzip extract --verbose {input.graph} {input.fasta} {input.depth} {output.segments} {output.depth} {output.fasta}"


rule megahit_assemble:
    output:
        dir=directory("data/group/{group}/r.{stem}.megahit-k{ksize}.d"),
        fasta="data/group/{group}/r.{stem}.megahit-k{ksize}.fn",
    wildcard_constraints:
        ksize=integer_wc,
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        r1=lambda w: ",".join(
            [
                f"data/reads/{mgen}/r1.{w.stem}.fq.gz"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
        r2=lambda w: ",".join(
            [
                f"data/reads/{mgen}/r2.{w.stem}.fq.gz"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
    threads: 36
    shell:
        """
        megahit -t {threads} --k-list {wildcards.ksize} -o {output.dir} -1 {params.r1} -2 {params.r2}
        ln {output.dir}/final.contigs.fa {output.fasta}
        """


rule quality_asses_assembly_against_one_ref:
    output:
        directory("{stem}.gquast-{genome}.d"),
    input:
        tigs="{stem}.fn",
        ref="data/genome/{genome}.fn",
    threads: 24
    params:
        min_tig_length=1000,
    conda:
        "conda/quast.yaml"
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_tig_length} -R {input.ref} --output-dir {output} {input.tigs}
        """


rule quality_asses_assembly_against_all_refs:
    output:
        dir=directory("{stem}.quast-{group}.d"),
        contig_to_genome="{stem}.quast-{group}.contig_to_genome.txt",
    input:
        tigs="{stem}.fn",
        refs=lambda w: [
            f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]
        ],
    threads: 12
    params:
        min_tig_length=500,
        refs=lambda w: ",".join(
            [f"data/genome/{genome}.fn" for genome in config["genome_group"][w.group]]
        ),
    conda:
        "conda/quast.yaml"
    shell:
        """
        metaquast.py --silent --fragmented --threads={threads} --min-contig {params.min_tig_length} -r {params.refs} --output-dir {output.dir} {input.tigs}
        cp {output.dir}/combined_reference/contigs_reports/alignments_*.tsv {output.contig_to_genome}
        # Reduce storage requirements:
        rm -rf {output.dir}/{{combined_reference/icarus_viewers,quast_corrected_input,runs_per_reference}}
        """


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

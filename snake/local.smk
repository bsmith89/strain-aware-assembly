rule link_local_data_directories:
    output:
        directory(config["local_data_dirs"]),
    input:
        [f"{config['local_data_root']}/{dir}" for dir in config["local_data_dirs"]],
    params:
        root=config["local_data_root"],
    shell:
        """
        for dir in {output}
        do
            ln -s "{params.root}/$dir"
        done
        """


rule link_secure_local_data_directories:
    output:
        directory(config["secure_local_data_dirs"]),
    input:
        [
            f"{config['secure_local_data_root']}/{dir}"
            for dir in config["secure_local_data_dirs"]
        ],
    params:
        root=config["secure_local_data_root"],
    shell:
        """
        for dir in {output}
        do
            ln -s "{params.root}/$dir"
        done
        """


rule link_xjin_reference_genome:
    output:
        "raw/genome/xjin/{stem}",
    input:
        "/pollard/data/microbial_genomes/fischbachBiohubStrains/Hybrid_Closed_Genomes/{stem}",
    shell:
        alias_recipe_norelative


rule link_processed_ucfmt_d97_data:
    output: "data/reads/DS0097_{stem}/{r}.proc.fq.gz"
    input: "/pollard/home/bsmith/Projects/strain-corr/data/reads/DS0097_{stem}/{r}.proc.fq.gz"
    shell:
        alias_recipe_norelative

ruleorder: link_processed_ucfmt_d97_data > alias_cleaned_reads

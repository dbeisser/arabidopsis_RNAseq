# build and run kraken for classification

rule kraken_db_maker:
    output:
        "{results}/classification/krakendb_nt/seqid2taxid.map",
        "{results}/classification/krakendb_nt/taxonomy/names.dmp"
    conda:
        "../envs/classify.yaml"
    params:
        dir = "{results}/classification/krakendb_nt"
    shell:
        """
            kraken2-build --download-taxonomy --threads 18 --db {params.dir}
            kraken2-build --download-library nt --threads 18 --db {params.dir}
            kraken2-build --build --threads 18 --db {params.dir}
        """

rule kraken_assign_tax:
    input:
        seq = expand("{{results}}/qc/cleaned/{{sample}}_{read}.{format}.gz", read = config["reads"], format = config["format"]) if config["format"] == "fastq" else expand("{{results}}/qc/cleaned/{{sample}}.{format}.gz", format = config["format"]),
        database = "{results}/classification/krakendb_nt/seqid2taxid.map"
    output:
        "{results}/classification/{sample}/retax_combined.kout"
    params:
        db = "{results}/classification/krakendb_nt"
    threads:
        config["threads"]
    conda:
        "../envs/classify.yaml"
    shell: "kraken2 --output {output} --threads {threads} --db {params.db} {input.seq}"

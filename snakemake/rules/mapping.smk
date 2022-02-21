from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# data for arabidopsis: https://plants.ensembl.org/info/data/ftp/index.html

rule download_transcriptome:
    input:
        FTP.remote(expand("{transcriptome}", transcriptome = config["transcriptome"]))
    output:
        "{results}/transcriptome/trans.fa.gz"
    shell: "curl {input} -o {output}"

rule index_transcriptome:
    input:
        "{results}/transcriptome/trans.fa.gz"
    output:
        "{results}/transcriptome/trans_index/versionInfo.json"
    params:
        out = "{results}/transcriptome/trans_index"
    conda:
        "../envs/mapping.yaml"
    shell: "salmon index -k 29 -t {input} -i {params.out}"

# for RNA-Seq data
# rule quantify:
#     input: index = "{results}/transcriptome/trans_index/versionInfo.json",
#         reads = expand("{{results}}/qc/cleaned/{{sample}}_{read}.{format}.gz", read = config["reads"], format = config["format"]) if config["format"] == "fastq" else expand("{{results}}/qc/cleaned/{{sample}}.{format}.gz", format = config["format"])
#     output: "{results}/quant/{sample}/quant.sf"
#     conda:
#         "../envs/mapping.yaml"
#     params:
#         index = "{results}/transcriptome/trans_index",
#         out = "{results}/quant/{sample}"
#     shell: "salmon quant -i {params.index} -l A -r {input.reads} -p 20 --validateMappings -o {params.out} --gcBias"

# for quantseq data according to: https://www.nature.com/articles/s41598-019-55434-x#Sec10
# We ran Salmon with the following default settings for the RNA-seq data: quant -i transcripts_index -l ISR -1 *R1*fastq* -2 *R2*fastq* -o $output -p 6–posBias–gcBias–seqBias –writeUnmappedNames. For the QuantSeq data which would only capture the 3′ end of the transcript we changed the settings to: quant -i transcripts_index -l SR -r *R1_trimmed.fq* -o $output -p 6 –writeUnmappedNames. The Bioconductor package tximport22 was used to convert transcripts to gene counts which could be used in downstream differential expression analysis.
# SR: wrong -> SF

rule quantify:
    input:
        index = "{results}/transcriptome/trans_index/versionInfo.json",
        reads = expand("{{results}}/qc/cleaned/{{sample}}_{read}.{format}.gz", read = config["reads"], format = config["format"]) if config["format"] == "fastq" else expand("{{results}}/qc/cleaned/{{sample}}.{format}.gz", format = config["format"])
    output:
        "{results}/quant/{sample}/quant.sf"
    conda:
        "../envs/mapping.yaml"
    params:
        index = "{results}/transcriptome/trans_index",
        out = "{results}/quant/{sample}"
    shell: "salmon quant -i {params.index} -l SF -r {input.reads} -p 20 -o {params.out}"

rule mapping_rate:
    input:
        "{results}/quant/{sample}/logs/salmon_quant.log"
    output:
        "{results}/quant/{sample}/logs/mapping_rate.csv"
    shell: "echo {wildcards.sample} $(grep 'Mapping rate =' {input}) > {output}"

rule comb_mapping_rate:
    input:
        expand("{{results}}/quant/{sample}/logs/mapping_rate.csv", sample=config["samples"])
    output:
        "{results}/quant/mapping_rates.csv"
    shell: "cat {input} > {output}"

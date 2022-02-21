# import modules
import os
import glob

# functions
def get_fastq(wildcards):
    files = glob.glob(os.path.join(config['raw'], wildcards.sample+'_*'+wildcards.read+'*.f*q.gz'))
    return sorted(files)

# merging and fastqc only for fastq files
if (config["format"] == "fastq"):
    # merge several lanes
    # not neccessary in this case, but also doesn't hurt
    rule merge:
        input: get_fastq
        output:
            merged = "{results}/merged_fastq/{sample}_{read}.fastq.gz",
            unzipped = temp("{results}/merged_fastq/{sample}_{read}.fastq")
        wildcard_constraints:
            read="R[1,2]"
        resources: io=1
        run:
            infiles = " ".join(input)
            print(infiles)
            shell("zcat {infiles} > {output.unzipped}")
            shell("gzip -c {output.unzipped} > {output.merged}")

    # fastqc analysis
    rule fastqc:
        input:
            "{results}/merged_fastq/{sample}_{read}.fastq.gz"
        output:
            zipped = temp("{results}/qc/fastqc/{sample}_{read}_fastqc.zip"),
            html = "{results}/qc/fastqc/{sample}_{read}_fastqc.html"
        conda:
            "../envs/preprocessing.yaml"
        params:
            outdir = "{results}/qc/fastqc"
        shell:
            """
                fastqc {input} --outdir={params.outdir}
            """

    # multiqc, summarizes fastqc reports
    rule multiqc:
        input:
            expand("{{results}}/qc/fastqc/{sample}_{read}_fastqc.zip", sample=config["samples"], read = config["reads"])
        output:
            "{results}/qc/multiqc_report.html"
        params:
            indir = "{results}/qc/fastqc",
            outdir = "{results}/qc",
            filename = "multiqc_report"
        conda:
            "../envs/preprocessing.yaml"
        shell:
            """
                multiqc {params.indir} -o {params.outdir} -n {params.filename}
            """

# cut adapters and low quality reads
rule cutadapt:
    input:
        expand("{{results}}/merged_fastq/{{sample}}_{read}.{format}.gz", read = config["reads"], format=config["format"]) if config["format"] == "fastq" else expand("{raw}/{{sample}}.{format}", raw=config["raw"], format=config["format"])
    output:
        expand("{{results}}/qc/cleaned/{{sample}}_{read}.{format}.gz", read = config["reads"], format = config["format"]) if config["format"] == "fastq" else expand("{{results}}/qc/cleaned/{{sample}}.{format}.gz", format = config["format"])
    log:
        "{results}/qc/cleaned/{sample}.log"
    params:
        adapter = config["adapter"],
        format = config["format"]
    threads: config["threads"]
    conda:
        "../envs/preprocessing.yaml"
    script:
        "../scripts/cutadapt.py"

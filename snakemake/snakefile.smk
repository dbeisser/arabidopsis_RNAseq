# Task: Analyse A. thaliana transcriptome samples
# Author: Daniela Beisser

configfile: "config.yaml"
report: "report/workflow.rst"

# rule all #####################################################################

rule all:
    input:
        expand("{results}/quant/{sample}/quant.sf", results=config["results"], sample=config["samples"]),
        expand("{results}/qc/multiqc_report.html", results=config["results"]),
        expand("{results}/quant/mapping_rates.csv", results=config["results"]),
        expand("{results}/classification/{sample}/retax_combined.kout", results=config["results"], sample=config["samples"])

# rules ########################################################################

# pre-processing
include: "rules/preprocessing.smk"

# quantification
include: "rules/mapping.smk"

# classify
include: "rules/classify.smk"

# status messages ##############################################################

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")
    shell("cat {log}")

################################################################################

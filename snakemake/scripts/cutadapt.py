import subprocess

logfile = open(str(snakemake.log), "w")
if len(snakemake.input) == 2:
    if snakemake.params["format"] == "fastq":
        subprocess.call(["cutadapt", "-a", str(snakemake.params.adapter), "-q 20", "-j", str(snakemake.threads), "-m 30", "-o", snakemake.output[0], "-p",
                                 snakemake.output[1], snakemake.input[0], snakemake.input[1]], stdout=logfile)
    else:
        subprocess.call(["cutadapt", "-a", str(snakemake.params.adapter), "-j", str(snakemake.threads), "-m 30", "-o", snakemake.output[0], "-p",
                                 snakemake.output[1], snakemake.input[0], snakemake.input[1]], stdout=logfile)
else:
    if snakemake.params["format"] == "fastq":
        subprocess.call(["cutadapt", "-a", str(snakemake.params.adapter), "-q 20", "-j", str(snakemake.threads), "-m 30", "-o", snakemake.output[0], snakemake.input[0]], stdout=logfile)
    else:
        subprocess.call(["cutadapt", "-a", str(snakemake.params.adapter), "-m 30", "-j", str(snakemake.threads), "-o", snakemake.output[0], snakemake.input[0]], stdout=logfile)
logfile.close()

# Load the configuration file
configfile: "config/config.yaml"

# Access the variables
RAW_DIR = config["directories"]["raw"]

# Sample accessions from GSE52778
# Untreated: SRR1039508, SRR1039512, SRR1039516, SRR1039520
# Dex-treated: SRR1039509, SRR1039513, SRR1039517, SRR1039521
SAMPLES = ['SRR1039508', 'SRR1039512', 'SRR1039516', 'SRR1039520', 'SRR1039509', 'SRR1039513', 'SRR1039517', 'SRR1039521']

rule all:
    input:
        expand(f"{RAW_DIR}/{{sample}}_1.fastq.gz", sample=SAMPLES)

include: "rules/01_download.smk"
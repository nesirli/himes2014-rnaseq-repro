# Load the configuration file
configfile: "config/config.yaml"

# Access the variables
RAW_DIR = config["directories"]["raw"]
REF_DIR = config["directories"]["reference"]

# Sample accessions from GSE52778
# Untreated: SRR1039508, SRR1039512, SRR1039516, SRR1039520
# Dex-treated: SRR1039509, SRR1039513, SRR1039517, SRR1039521
SAMPLES = ['SRR1039508', 'SRR1039512', 'SRR1039516', 'SRR1039520', 'SRR1039509', 'SRR1039513', 'SRR1039517', 'SRR1039521']

rule all:
    input:
        expand(f"{RAW_DIR}/{{sample}}_1.fastq.gz", sample=SAMPLES),
        f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf",

include: "rules/01_download.smk"
# include: "rules/02_index.smk"
# include: "rules/03_qc_trim.smk"
# include: "rules/04_align.smk"
# include: "rules/05_counts.smk"
# include: "rules/06_multiqc.smk"
# include: "rules/07_deseq2.smk"
# include: "rules/08_figures.smk"
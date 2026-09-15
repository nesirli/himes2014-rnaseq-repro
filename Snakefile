# Load the configuration file
configfile: "config/config.yaml"


# Access the variables
RAW_DIR = config["directories"]["raw"]
REF_DIR = config["directories"]["reference"]
QC_DIR = config["directories"]["qc"]
TRIM_DIR = config["directories"]["trim"]
ALIGN_DIR = config["directories"]["align"]

# Sample accessions from GSE52778
# Untreated: SRR1039508, SRR1039512, SRR1039516, SRR1039520
# Dex-treated: SRR1039509, SRR1039513, SRR1039517, SRR1039521
SAMPLES = [
    "SRR1039508",
    "SRR1039512",
    "SRR1039516",
    "SRR1039520",
    "SRR1039509",
    "SRR1039513",
    "SRR1039517",
    "SRR1039521",
]


rule all:
    input:
        expand(f"{RAW_DIR}/{{sample}}_1.fastq.gz", sample=SAMPLES),
        f"{RAW_DIR}/data_stats.txt",
        f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf",
        f"{REF_DIR}/reference_integrity.txt",
        f"{REF_DIR}/gene_id2name.tsv",
        f"{QC_DIR}/post_trim/multiqc_report.html",
        expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai", sample=SAMPLES),
        "results/counts/gene_counts.txt"


include: "rules/01_download.smk"
include: "rules/02_qc_trim.smk"
include: "rules/03_index_align.smk"
include: "rules/04_counts.smk"
# include: "rules/05_deseq2.smk"
# include: "rules/06_figures.smk"

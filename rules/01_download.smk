rule download_samples:
    output:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
    log:
        "logs/download/{sample}.log"
    conda:
        "../envs/01_download.yaml"
    threads:
        config["params"]["download-threads"]
    shell:
        """
        # Download raw data files
        fasterq-dump {wildcards.sample} --split-files --threads {threads} --outdir {RAW_DIR} > {log} 2>&1
        gzip {RAW_DIR}/{wildcards.sample}_1.fastq >> {log} 2>&1
        gzip {RAW_DIR}/{wildcards.sample}_2.fastq >> {log} 2>&1
        """

rule download_reference:
    output:
        genome=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf=f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf"
    log:
        "logs/download/reference.log"
    conda:
        "../envs/01_download.yaml"
    threads:
        config["params"]["download-threads"]
    shell:
        """
        wget -P {REF_DIR} https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > {log} 2>&1
        gunzip {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz >> {log} 2>&1

        wget -P {REF_DIR} https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz >> {log} 2>&1
        gunzip {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz >> {log} 2>&1
        """
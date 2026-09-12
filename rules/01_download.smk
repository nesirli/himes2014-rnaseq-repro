rule download_samples:
    output:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        stats=f"{RAW_DIR}/{{sample}}_stats.txt"
    log:
        "logs/download/{sample}.log"
    conda:
        "../envs/01_download.yaml"
    threads:
        config["params"]["download-threads"]
    shell:
        """
        set -euo pipefail

        # Prefetch the SRA object first
        prefetch {wildcards.sample} --output-directory {RAW_DIR}/sra_cache > {log} 2>&1

        # Validate checksums against NCBI's stored values before conversion
        vdb-validate {RAW_DIR}/sra_cache/{wildcards.sample}/{wildcards.sample}.sra >> {log} 2>&1

        # Convert to FASTQ, --split-3 keeps _1/_2 counts consistent
        fasterq-dump {RAW_DIR}/sra_cache/{wildcards.sample}/{wildcards.sample}.sra \
            --split-3 --threads {threads} --outdir {RAW_DIR} >> {log} 2>&1

        # Compress paired reads
        gzip {RAW_DIR}/{wildcards.sample}_1.fastq >> {log} 2>&1
        gzip {RAW_DIR}/{wildcards.sample}_2.fastq >> {log} 2>&1

        # Handle orphan/unpaired reads if --split-3 produced any
        if [ -f {RAW_DIR}/{wildcards.sample}.fastq ]; then
            gzip {RAW_DIR}/{wildcards.sample}.fastq >> {log} 2>&1
        fi

        # Per-sample stats, one header, no shared-file race condition
        seqkit stats {output.r1} {output.r2} > {output.stats}
        """

rule merge_stats:
    input:
        expand(f"{RAW_DIR}/{{sample}}_stats.txt", sample=SAMPLES)
    output:
        f"{RAW_DIR}/data_stats.txt"
    shell:
        """
        set -euo pipefail
        head -n 1 {input[0]} > {output}
        for f in {input}; do
            tail -n +2 "$f" >> {output}
        done
        """

rule download_reference:
    output:
        genome=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf=f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf",
        integrity=f"{REF_DIR}/reference_integrity.txt"
    log:
        "logs/download/reference.log"
    conda:
        "../envs/01_download.yaml"
    threads:
        config["params"]["download-threads"]
    shell:
        """
        set -euo pipefail

        # Download the archives and Ensembl's checksum manifests separately
        # (both CHECKSUMS endpoints are literally named "CHECKSUMS", so keep them apart)
        wget -O {REF_DIR}/CHECKSUMS_dna https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/CHECKSUMS >> {log} 2>&1
        wget -O {REF_DIR}/CHECKSUMS_gtf https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/CHECKSUMS >> {log} 2>&1

        wget -O {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz >> {log} 2>&1
        wget -O {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz >> {log} 2>&1

        # Record local BSD checksums alongside the expected Ensembl values
        sum {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > {output.integrity}
        grep "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" {REF_DIR}/CHECKSUMS_dna >> {output.integrity}
        sum {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz >> {output.integrity}
        grep "Homo_sapiens.GRCh38.110.gtf.gz" {REF_DIR}/CHECKSUMS_gtf >> {output.integrity}

        gunzip -f {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz >> {log} 2>&1
        gunzip -f {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz >> {log} 2>&1
        """
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
    resources:
        mem_mb=2000
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
    resources:
        mem_mb=1000
    shell:
        """
        set -euo pipefail
        head -n 1 {input[0]} > {output}
        for f in {input}; do
            tail -n +2 "$f" >> {output}
        done
        """

rule gene_id2name:
    input:
        gtf=f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf"
    output:
        f"{REF_DIR}/gene_id2name.tsv"
    log:
        "logs/download/gene_id2name.log"
    resources:
        mem_mb=1000
    shell:
        """
        awk 'BEGIN {{ FS = OFS = "\t" }}
             $3 == "gene" {{
                 id = ""; nm = ""
                 if (match($9, /gene_id "[^"]+"/))   id = substr($9, RSTART + 9,  RLENGTH - 10)
                 if (match($9, /gene_name "[^"]+"/)) nm = substr($9, RSTART + 11, RLENGTH - 12)
                 if (id != "") print id, nm
             }}' {input.gtf} > {output} 2> {log}
        """

rule download_reference:
    output:
        gtf=f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf",
        integrity=f"{REF_DIR}/reference_integrity.txt"
    log:
        "logs/download/reference.log"
    conda:
        "../envs/01_download.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        set -euo pipefail

        # Download the archive and Ensembl's checksum manifest
        wget -O {REF_DIR}/CHECKSUMS_gtf https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/CHECKSUMS >> {log} 2>&1
        wget -O {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz >> {log} 2>&1

        # Record the local BSD checksum alongside the expected Ensembl value
        sum {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz > {output.integrity}
        grep "Homo_sapiens.GRCh38.110.gtf.gz" {REF_DIR}/CHECKSUMS_gtf >> {output.integrity}

        gunzip -f {REF_DIR}/Homo_sapiens.GRCh38.110.gtf.gz >> {log} 2>&1
        """

rule download_index:
    output:
        multiext(f"{HISAT2_INDEX_DIR}/genome",
                ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")
    log:
        "logs/download/download_index.log"
    conda:
        "../envs/01_download.yaml"
    resources:
        mem_mb=2000
    params:
        url=config["params"]["hisat2-index-url"],
        index_dir=HISAT2_INDEX_DIR
    shell:
        """
        set -euo pipefail
        mkdir -p {params.index_dir}
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        curl -L --fail --retry 5 --retry-delay 5 --no-progress-meter \
            -o "$tmp/index.tar.gz" {params.url} >> {log} 2>&1
        tar -xzf "$tmp/index.tar.gz" -C "$tmp" >> {log} 2>&1
        n=$(find "$tmp" -name '*.ht2' | wc -l | tr -d ' ')
        if [ "$n" -ne 8 ]; then
            echo "expected 8 .ht2 files, found $n" >&2
            exit 1
        fi
        find "$tmp" -name '*.ht2' | while read -r f; do
            mv "$f" "{params.index_dir}/genome.$(basename "$f" | cut -d. -f2-)"
        done
        """

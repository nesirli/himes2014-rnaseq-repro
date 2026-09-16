HISAT2_INDEX_DIR = f"{REF_DIR}/hisat2_index/prebuilt"

rule download_index:
    output:
        multiext(f"{HISAT2_INDEX_DIR}/genome",
                ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")
    log:
        "logs/index_align/download_index.log"
    conda:
        "../envs/03_index_align.yaml"
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

rule hisat2_align:
    input:
        r1=f"{TRIM_DIR}/{{sample}}_1_trim.fastq.gz",
        r2=f"{TRIM_DIR}/{{sample}}_2_trim.fastq.gz",
        index=multiext(f"{HISAT2_INDEX_DIR}/genome",
                ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")
    output:
        bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam",
        bai=f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai"
    log:
        "logs/index_align/{sample}.log"
    conda:
        "../envs/03_index_align.yaml"
    threads:
        config["params"]["align-threads"]
    resources:
        mem_mb=16000
    params:
        index=f"{HISAT2_INDEX_DIR}/genome"
    shell:
        """
        mkdir -p {ALIGN_DIR}
        hisat2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} 2>> {log} \
            | samtools sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """
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
        "logs/align/{sample}.log"
    conda:
        "../envs/03_align.yaml"
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

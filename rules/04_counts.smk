rule feature_counts:
    input:
        bam=expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
        gtf=f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf"
    output:
        r="results/counts/gene_counts.txt"
    log: 
        "logs/counts/counts.log"
    conda:
        "../envs/04_counts.yaml"
    threads:
        config["params"]["count-threads"]
    resources:
        mem_mb=8000
    shell:
        """
        featureCounts \
            -a {input.gtf} \
            -o {output.r} \
            -T {threads} \
            -p --countReadPairs \
            -s 2 \
            {input.bam} > {log} 2>&1
        """

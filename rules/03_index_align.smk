rule hisat2_index:
    input:
        genome=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf=f"{REF_DIR}/Homo_sapiens.GRCh38.110.gtf"
    output:
        multiext(f"{REF_DIR}/hisat2_index/genome",
                ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")
    log: 
        "logs/index_align/hisat2_index.log"
    conda:
        "../envs/03_index_align.yaml"
    threads:
        config["params"]["index-threads"]
    resources:
        mem_mb=200000
    params:
        index_dir=f"{REF_DIR}/hisat2_index"
    shell:
        """
        mkdir -p {params.index_dir}
        hisat2_extract_splice_sites.py {input.gtf} > {params.index_dir}/genome.ss >> {log} 2>&1
        hisat2_extract_exons.py {input.gtf} > {params.index_dir}/genome.exon >> {log} 2>&1
        hisat2-build -p {threads} \
            --ss {params.index_dir}/genome.ss \
            --exon {params.index_dir}/genome.exon \
            {input.genome} {params.index_dir}/genome >> {log} 2>&1
        """

rule hisat2_align:
    input:
        r1=f"{TRIM_DIR}/{{sample}}_1_trim.fastq.gz",
        r2=f"{TRIM_DIR}/{{sample}}_2_trim.fastq.gz",
        index=multiext(f"{REF_DIR}/hisat2_index/genome",
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
        index=f"{REF_DIR}/hisat2_index/genome"
    shell:
        """
        mkdir -p {ALIGN_DIR}
        hisat2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} 2>> {log} \
            | samtools sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """
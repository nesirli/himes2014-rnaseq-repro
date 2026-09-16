rule pre_trim_fastqc:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
    output:
        html1=f"{QC_DIR}/pre_trim/{{sample}}_1_fastqc.html",
        html2=f"{QC_DIR}/pre_trim/{{sample}}_2_fastqc.html",
    log:
        "logs/qc/pre_trim/{sample}.log"
    conda:
        "../envs/02_qc_trim.yaml"
    threads:
        config["params"]["qc-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        mkdir -p {QC_DIR}/pre_trim
        fastqc {input.r1} {input.r2} -o {QC_DIR}/pre_trim --threads {threads} > {log} 2>&1
        """

rule pre_trim_multi_qc:
    input:
        expand(f"{QC_DIR}/pre_trim/{{sample}}_1_fastqc.html", sample=SAMPLES),
        expand(f"{QC_DIR}/pre_trim/{{sample}}_2_fastqc.html", sample=SAMPLES),
    output:
        f"{QC_DIR}/pre_trim/multiqc_report.html"
    log:
        "logs/qc/pre_trim/multiqc.log"
    conda:
        "../envs/02_qc_trim.yaml"
    threads:
        config["params"]["qc-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        multiqc {QC_DIR}/pre_trim -o {QC_DIR}/pre_trim > {log} 2>&1
        """

rule trim:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        report=f"{QC_DIR}/pre_trim/multiqc_report.html"
    output:
        t1=f"{TRIM_DIR}/{{sample}}_1_trim.fastq.gz",
        t2=f"{TRIM_DIR}/{{sample}}_2_trim.fastq.gz",
        html=f"{TRIM_DIR}/{{sample}}_fastp.html",
        json=f"{TRIM_DIR}/{{sample}}_fastp.json",
    log:
        "logs/trim/{sample}.log"
    conda:
        "../envs/02_qc_trim.yaml"
    threads:
        config["params"]["trim-threads"]
    resources:
        mem_mb=4000
    params:
        min_qual=config["params"]["phred-cutoff"]
    shell:
        """
        mkdir -p {TRIM_DIR}
        fastp -w {threads} \
            -i {input.r1} -I {input.r2} \
            -o {output.t1} -O {output.t2} \
            -q {params.min_qual} \
            -h {output.html} -j {output.json} \
            > {log} 2>&1
        """

rule post_trim_fastqc:
    input:
        t1=f"{TRIM_DIR}/{{sample}}_1_trim.fastq.gz",
        t2=f"{TRIM_DIR}/{{sample}}_2_trim.fastq.gz",
    output:
        html1=f"{QC_DIR}/post_trim/{{sample}}_1_trim_fastqc.html",
        html2=f"{QC_DIR}/post_trim/{{sample}}_2_trim_fastqc.html",
    log:
        "logs/qc/post_trim/{sample}.log"
    conda:
        "../envs/02_qc_trim.yaml"
    threads:
        config["params"]["qc-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        mkdir -p {QC_DIR}/post_trim
        fastqc {input.t1} {input.t2} -o {QC_DIR}/post_trim --threads {threads} > {log} 2>&1
        """

rule post_trim_multi_qc:
    input:
        expand(f"{QC_DIR}/post_trim/{{sample}}_1_trim_fastqc.html", sample=SAMPLES),
        expand(f"{QC_DIR}/post_trim/{{sample}}_2_trim_fastqc.html", sample=SAMPLES),
    output:
        f"{QC_DIR}/post_trim/multiqc_report.html"
    log:
        "logs/qc/post_trim/multiqc.log"
    conda:
        "../envs/02_qc_trim.yaml"
    threads:
        config["params"]["qc-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        multiqc {QC_DIR}/post_trim -o {QC_DIR}/post_trim > {log} 2>&1
        """

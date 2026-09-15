rule deseq2:
    input:
        counts="results/counts/gene_counts.txt",
        gene_map=f"{REF_DIR}/gene_id2name.tsv"
    output:
        all="results/de/all_results.csv",
        sig="results/de/significant_genes.csv"
    log:
        "logs/deseq2/deseq2.log"
    conda:
        "../envs/05_deseq2.yaml"
    script:
        "../scripts/01_deseq2.R"

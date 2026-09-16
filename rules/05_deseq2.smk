rule deseq2:
    input:
        counts="results/counts/gene_counts.txt",
        gene_map=f"{REF_DIR}/gene_id2name.tsv"
    output:
        all_results="results/de/all_results.csv",
        sig="results/de/significant_genes.csv",
        volcano="results/figures/volcano.png",
        pca="results/figures/pca.png",
        heatmap="results/figures/heatmap_top30.png"
    log:
        "logs/deseq2/deseq2.log"
    resources:
        mem_mb=8000
    conda:
        "../envs/05_deseq2.yaml"
    script:
        "../scripts/01_deseq2.R"

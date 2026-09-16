# Himes 2014 RNA-seq reproduction

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A57-039475)](https://snakemake.readthedocs.io/)
[![Conda](https://img.shields.io/badge/Conda-environments-44A833?logo=anaconda&logoColor=white)](https://docs.conda.io/)
[![HISAT2](https://img.shields.io/badge/HISAT2-aligner-1f6feb)](https://daehwankimlab.github.io/hisat2/)
[![featureCounts](https://img.shields.io/badge/featureCounts-subread-555555)](https://subread.sourceforge.net/)
[![DESeq2](https://img.shields.io/badge/DESeq2-Bioconductor-276DC3)](https://bioconductor.org/packages/DESeq2/)
[![R](https://img.shields.io/badge/R-Bioconductor-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)

Reproduction of Himes et al. (2014), *RNA-Seq Transcriptome Profiling Identifies
CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in
Airway Smooth Muscle Cells*, using a Snakemake pipeline that downloads the raw
reads, performs QC/trimming, aligns them with HISAT2, quantifies with
featureCounts, and calls differential expression (dexamethasone vs. untreated)
with DESeq2.

## Data

- **GEO accession:** [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)
- **Samples (paired-end):**
  - Untreated: `SRR1039508`, `SRR1039512`, `SRR1039516`, `SRR1039520`
  - Dexamethasone-treated: `SRR1039509`, `SRR1039513`, `SRR1039517`, `SRR1039521`
- **Reference:** Ensembl *Homo sapiens* GRCh38, release 110 GTF
  (checksum-verified on download; see `data/reference/reference_integrity.txt`)
- **Aligner index:** prebuilt HISAT2 GRCh38 `genome_tran` index (transcript-aware
  graph index) from the [AWS Public Datasets](https://registry.opendata.aws/jhu-indexes/),
  downloaded rather than built locally

## Pipeline

1. **Download** — `prefetch`/`fasterq-dump` retrieve the SRA runs and
   `vdb-validate` checks integrity; the Ensembl GTF is verified against Ensembl's
   `CHECKSUMS` manifest; the prebuilt HISAT2 `genome_tran` index is fetched from
   AWS.
2. **QC & trimming** — FastQC before and after `fastp` trimming, summarized with
   MultiQC.
3. **Alignment** — HISAT2 against the prebuilt `genome_tran` graph index +
   `samtools`.
4. **Quantification** — `featureCounts` (paired-end, reverse-stranded).
5. **Differential expression** — DESeq2 (`~condition`), with figures.

## Tools and versions

| Tool | Version | Stage |
| --- | --- | --- |
| Snakemake | 9.26.1 | workflow |
| sra-tools | 3.4.1 | download |
| wget | 1.25.0 | download |
| seqkit | 2.13.0 | download |
| FastQC | 0.12.1 | QC |
| fastp | 0.24.3 | trimming |
| MultiQC | 1.35 | QC |
| HISAT2 | 2.2.3 | alignment |
| SAMtools | 1.24 | alignment |
| featureCounts (Subread) | 2.1.1 | quantification |
| R | 4.5.3 | differential expression |
| DESeq2 | 1.50.2 | differential expression |
| ggplot2 | 4.0.3 | figures |
| pheatmap | 1.0.13 | figures |

The exact versions are pinned in the `envs/` files. The reference is Ensembl
GRCh38 release 110, and the aligner uses the prebuilt HISAT2 `genome_tran` index.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) (>= 7)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)
- Per-stage, pinned environments are provided in `envs/` (sra-tools, FastQC,
  fastp, MultiQC, HISAT2, samtools, subread, R/Bioconductor/DESeq2, ggplot2,
  pheatmap).

## Usage

Run the full pipeline from the repository root:

```bash
snakemake --use-conda --cores 8
```

Thread counts and the Phred cutoff can be adjusted in `config/config.yaml`.
Conda environments are created automatically on first run.

`--cores` is the **total** CPU budget, not a per-job cap, and must be at least
the largest per-rule value in `config/config.yaml` (currently 8 for
`align-threads`). Set it to your machine's core count to run independent samples
in parallel — e.g. `--cores 16` runs two 8-thread alignments at once.

Rules also declare `resources: mem_mb`, which Snakemake enforces if you pass
`--resources mem_mb=<N>`.

The HISAT2 index is **downloaded** (~4.1 GB tarball) instead of built locally, so
no large-memory machine is required. The index URL is configurable via
`params.hisat2-index-url` in `config/config.yaml`.

## Outputs

| Output | Path |
| --- | --- |
| Gene counts | `results/counts/gene_counts.txt` |
| All DE results | `results/de/all_results.csv` |
| Significant genes | `results/de/significant_genes.csv` |
| Volcano plot | `results/figures/volcano.png` |
| PCA plot | `results/figures/pca.png` |
| Top-30 heatmap | `results/figures/heatmap_top30.png` |

## Results

The run reproduces the main results of the paper.

- **Alignment:** 97.7–98.9% overall alignment rate for all 8 samples. All BAM
  files pass `samtools quickcheck`.
- **Counting:** 7.6–15.4 million reads were assigned to genes per sample. The
  unassigned categories look normal.
- **Differential expression:** 835 significant genes — 463 up and 372 down.
- **Validation:** the five known dexamethasone-responsive genes are all
  significant (padj < 0.05, log2FC > 1): `CRISPLD2`, `DUSP1`, `KLF15`, `PER1`
  and `TSC22D3`.
- **Figures:** the volcano plot, PCA plot and top-30 heatmap are in
  `results/figures/`. `results/de/significant_genes.csv` has 835 rows.

## Data and reference integrity

The pipeline writes simple text files so you can check the inputs:

- `data/raw/{sample}_stats.txt` — read count and read length for each downloaded
  FASTQ file. `data/raw/data_stats.txt` joins all samples in one table.
- `data/reference/reference_integrity.txt` — two lines for the Ensembl GTF: the
  local checksum and the checksum from Ensembl's `CHECKSUMS` file. When the two
  values match, the reference is complete and correct.
- `results/counts/gene_counts.txt.summary` — for each sample, how many reads
  `featureCounts` assigned to genes and how many it did not assign.

The pipeline also runs `vdb-validate` on every SRA file. This compares each file
with NCBI's stored checksums before converting it to FASTQ.

## Repository structure

```
Snakefile            # Entry point; wires the stages together
config/config.yaml   # Paths and tunable parameters
rules/               # One .smk file per pipeline stage
scripts/             # R analysis script (DESeq2 + figures)
envs/                # Pinned conda environments
data/                # Raw reads and reference (downloaded, git-ignored)
results/             # Pipeline outputs (git-ignored)
paper/               # Original article and supplementary data
```

## Citation

Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, et al. RNA-Seq
Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene
that Modulates Cytokine Function in Airway Smooth Muscle Cells. *PLoS ONE*
9(6): e99625 (2014). doi:[10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625)

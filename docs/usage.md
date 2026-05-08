# eDNA Metabarcoding Pipeline — Usage Guide

## Overview

Multi-marker eDNA pipeline supporting 16S, 18S, ITS, CO1, and 12S amplicon
sequencing from raw paired-end reads through to ecological analysis outputs.

### Pipeline steps

```
Raw reads (FASTQ)
    │
    ├── FastQC (raw QC)
    │
    ├── Cutadapt (primer trimming, per-marker primers)
    │
    ├── DADA2 error model learning (per sequencing run)
    │
    ├── DADA2 denoising → ASV table + sequences (per sample)
    │
    ├── VSEARCH chimera detection (de novo)
    │
    ├── Merge ASV tables across samples + final bimera removal
    │
    ├── Taxonomy assignment
    │   ├── 16S  → SILVA 138.1 (DADA2 naive Bayesian)
    │   ├── 18S  → PR2 v5.0 (DADA2 naive Bayesian)
    │   ├── ITS  → UNITE v10 (DADA2 naive Bayesian)
    │   ├── CO1  → MIDORI2 GB262 (DADA2 naive Bayesian)
    │   └── 12S  → MIDORI2 GB262 (DADA2 naive Bayesian)
    │
    ├── MultiQC report
    │
    ├── Basic ecological analysis (phyloseq + vegan)
    │   ├── Alpha diversity (Observed, Chao1, Shannon, Simpson)
    │   ├── Beta diversity (Bray-Curtis, Jaccard)
    │   ├── PCoA, NMDS ordination
    │   └── Taxonomic composition barplots (Phylum → Genus)
    │
    └── Full ecological analysis suite (9 modules, runs after basic)
        ├── 01 Alpha diversity — richness, evenness, rarefaction (iNEXT),
        │       rank-abundance, occupancy-abundance, Kruskal-Wallis + Dunn's
        ├── 02 Beta diversity — 6 distance metrics, PERMANOVA (single + multi-factor),
        │       pairwise PERMANOVA (BH-adjusted), PERMDISP, ANOSIM, MRPP
        ├── 03 Ordination — PCoA (Bray-Curtis + Aitchison), NMDS + Shepard diagram,
        │       PCA-CLR biplot, RDA with forward selection, db-RDA, CCA
        ├── 04 Differential abundance — DESeq2 (all pairwise) + ALDEx2,
        │       volcano plots, MA plots, significant-taxa heatmap
        ├── 05 Co-occurrence network — Spearman on CLR, significance threshold,
        │       Louvain modularity, hub/keystone taxa, ggraph visualisation
        ├── 06 Indicator species — IndVal (multipatt), SIMPER (all pairs),
        │       core microbiome (multiple thresholds), per-group taxon stats
        ├── 07 Environmental drivers — envfit vectors, Mantel test (per variable),
        │       variance partitioning, Procrustes (community vs env PCA),
        │       alpha diversity × env Spearman correlations
        ├── Cross-marker — Procrustes between markers, genus sharing matrix,
        │       multi-marker combined PCoA
        └── HTML report — per-marker R Markdown report integrating all results
```

---

## Quick start

### 1. Set up databases

```bash
bash assets/download_databases.sh databases/
```

Update `nextflow.config` database paths to match your local files.

### 2. Prepare samplesheet

See `assets/samplesheet_example.csv`. Required columns:

| Column    | Description |
|-----------|-------------|
| `sample`  | Sample identifier (no spaces) |
| `fastq_1` | Absolute path to R1 FASTQ (gzipped) |
| `fastq_2` | Absolute path to R2 FASTQ (gzipped). Leave blank for single-end |
| `marker`  | One of: `16S`, `18S`, `ITS`, `CO1`, `12S` |
| `run`     | Sequencing run ID (used to group error model learning). Optional — defaults to `run1` |

The same sample can appear multiple times with different `marker` values (one row per marker per sample).

### 3. Run the pipeline

**Single marker (16S), Docker:**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --markers 16S \
    -profile docker
```

**Multiple markers, Singularity on SLURM:**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --markers 16S,ITS,CO1 \
    --metadata metadata.tsv \
    -profile slurm,singularity \
    -resume
```

**All five markers with custom resources:**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --markers 16S,18S,ITS,CO1,12S \
    --metadata metadata.tsv \
    --max_cpus 32 \
    --max_memory 256.GB \
    -profile singularity \
    -resume
```

**Skip ecology (just produce ASV tables + taxonomy):**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --markers 16S \
    --run_ecology false \
    -profile docker
```

---

## Parameters

### Main pipeline

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Path to samplesheet CSV |
| `--outdir` | `results` | Output directory |
| `--markers` | `16S` | Comma-separated marker list |
| `--metadata` | null | Sample metadata TSV (optional, enhances ecology) |
| `--run_ecology` | `true` | Run basic ecological analysis |
| `--run_full_ecology` | `true` | Run full 9-module ecological analysis suite |
| `--dada2_pool` | `false` | DADA2 pooling: `false`, `true`, or `pseudo` |
| `--dada2_chimera` | `consensus` | Chimera method: `consensus`, `pooled`, `per-sample` |
| `--max_cpus` | `32` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum walltime per process |

### Full ecology options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ecology_group_var` | null | Metadata column to use as primary group (e.g. `habitat`). Auto-detects first categorical column if null |
| `--ecology_min_prevalence` | `0.3` | Network: minimum fraction of samples an ASV must appear in |
| `--ecology_cor_cutoff` | `0.6` | Network: minimum Spearman \|r\| for a co-occurrence edge |

### Overriding primers

Default primers are defined in `nextflow.config` under `params.primers`. Override per-run:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --markers 16S \
    --primers.16S.fwd GTGYCAGCMGCCGCGGTAA \
    --primers.16S.rev GGACTACNVGGGTWTCTAAT \
    --primers.16S.trunc_len_f 240 \
    --primers.16S.trunc_len_r 180 \
    -profile docker
```

---

## Output structure

```
results/
├── pipeline_info/          # Execution report, timeline, trace, DAG
├── fastqc/
│   └── {marker}/           # Per-sample FastQC reports
├── trimmed/
│   └── {marker}/           # Cutadapt logs (trimmed reads not retained by default)
├── multiqc/                # Aggregated QC report
├── dada2/
│   └── {marker}/
│       ├── error_models/   # DADA2 error model RDS + plots
│       └── asv_tables/     # Per-sample ASV tables + sequences + read stats
├── chimera_check/
│   └── {marker}/           # VSEARCH chimera detection results
├── asv_tables/
│   └── {marker}/
│       ├── *.merged_asv_table.tsv    # Final ASV count table (ASVs × samples)
│       ├── *.merged_asv_table.rds    # Phyloseq-ready RDS
│       └── *.read_tracking.tsv       # Reads surviving each step
├── taxonomy/
│   └── {marker}/
│       └── *.taxonomy.tsv            # ASV taxonomy assignments + bootstrap
└── ecology/
    └── {marker}/
        ├── diversity/
        │   ├── alpha_diversity.tsv
        │   ├── alpha_diversity.pdf
        │   ├── bray_curtis_distance.tsv
        │   ├── jaccard_distance.tsv
        │   ├── permanova_bray.tsv
        │   ├── anosim_bray.txt
        │   └── phyloseq_object.rds
        ├── ordination/
        │   ├── pcoa_bray_curtis.pdf
        │   ├── nmds_bray_curtis.pdf
        │   ├── pcoa_eigenvalues.tsv
        │   └── rda_triplot.pdf
        └── composition/
            ├── Phylum_barplot.pdf
            ├── Class_barplot.pdf
            ├── Order_barplot.pdf
            ├── Family_barplot.pdf
            ├── Genus_barplot.pdf
            ├── *_relative_abundance.tsv
            └── top_asvs_heatmap.pdf
```

---

## Running the full ecology pipeline standalone

If the main pipeline has already completed, run just the ecological analysis
on the existing results without re-processing reads:

```bash
# All markers, from previous results directory
nextflow run ecology_pipeline.nf \
    --results_dir results/ \
    --markers 16S,18S,ITS,CO1,12S \
    --metadata metadata.tsv \
    --ecology_group_var habitat \
    --outdir full_ecology/ \
    -profile singularity \
    -resume

# Single marker, from direct file inputs
nextflow run ecology_pipeline.nf \
    --asv_table results/asv_tables/16S/16S.merged_asv_table.tsv \
    --taxonomy  results/taxonomy/16S/SAMPLE.taxonomy.tsv \
    --marker    16S \
    --metadata  metadata.tsv \
    --outdir    full_ecology/ \
    -profile docker

# Skip full ecology during main run (process reads only):
nextflow run main.nf \
    --input samplesheet.csv \
    --markers 16S \
    --run_ecology false \
    --run_full_ecology false \
    -profile singularity
```

### Full ecology output structure

```
full_ecology/ (or results/full_ecology/ when run inline)
├── {MARKER}/
│   ├── 00_report/
│   │   └── {MARKER}_ecological_report.html   ← Start here
│   ├── 01_alpha_diversity/
│   │   ├── alpha_diversity_metrics.tsv
│   │   ├── alpha_diversity_boxplots.pdf
│   │   ├── alpha_kruskal_wallis.tsv
│   │   ├── alpha_dunn_posthoc.tsv
│   │   ├── rarefaction_curves.pdf
│   │   ├── sample_coverage.pdf
│   │   ├── rank_abundance.pdf
│   │   └── occupancy_abundance.tsv + .pdf
│   ├── 02_beta_diversity/
│   │   ├── {method}_distance.tsv              ← 6 distance matrices
│   │   ├── bray_curtis_heatmap.pdf
│   │   ├── permanova_results.tsv
│   │   ├── pairwise_permanova.tsv
│   │   ├── betadisper_permdisp.tsv
│   │   ├── anosim_results.txt
│   │   ├── mrpp_results.txt
│   │   └── multivariate_tests_summary.tsv
│   ├── 03_ordination/
│   │   ├── pcoa_bray_curtis.pdf + .png
│   │   ├── pcoa_aitchison.pdf + .png
│   │   ├── nmds_bray_curtis.pdf + .png
│   │   ├── nmds_shepard.pdf
│   │   ├── pca_clr_biplot.pdf + .png
│   │   ├── rda_environmental.pdf + .png
│   │   ├── rda_forward_selection.txt
│   │   ├── rda_anova.txt
│   │   ├── dbrda_bray_curtis.pdf + .png
│   │   ├── cca_environmental.pdf
│   │   └── ordination_panel.pdf + .png       ← 4-panel summary
│   ├── 04_differential_abundance/
│   │   ├── deseq2_{X}_vs_{Y}.tsv             ← All pairwise
│   │   ├── volcano_deseq2_{X}_vs_{Y}.pdf + .png
│   │   ├── ma_plot_deseq2.pdf
│   │   ├── aldex2_results.tsv
│   │   ├── aldex2_volcano.pdf + .png
│   │   └── deseq2_significant_heatmap.pdf
│   ├── 05_co_occurrence_network/
│   │   ├── correlation_matrix.tsv
│   │   ├── network_edges.tsv
│   │   ├── network_statistics.tsv
│   │   ├── node_metrics.tsv
│   │   ├── hub_keystone_taxa.tsv
│   │   ├── modularity.tsv
│   │   ├── network_plot.pdf + .png
│   │   └── edge_type_pie.pdf
│   ├── 06_indicator_species/
│   │   ├── indval_summary.txt
│   │   ├── indval_significant.tsv
│   │   ├── indval_barplot.pdf + .png
│   │   ├── simper_summary.txt
│   │   ├── simper_{X}_vs_{Y}.tsv
│   │   ├── core_microbiome_matrix.tsv
│   │   ├── core_microbiome_plot.pdf + .png
│   │   └── per_group_taxon_stats.tsv
│   └── 07_env_drivers/
│       ├── envfit_summary.txt
│       ├── envfit_vectors.tsv
│       ├── nmds_envfit_vectors.pdf + .png
│       ├── mantel_tests.tsv
│       ├── mantel_results.pdf
│       ├── variance_partitioning.txt + .pdf + .tsv
│       ├── procrustes_analysis.txt
│       ├── procrustes_stats.tsv
│       ├── procrustes_plot.pdf + .png
│       └── alpha_env_correlations.tsv + heatmap.pdf
└── cross_marker/
    ├── crossmarker_alpha_diversity.tsv
    ├── alpha_diversity_by_marker.pdf + .png
    ├── procrustes_between_markers.tsv
    ├── marker_correlation_heatmap.pdf
    ├── genus_jaccard_between_markers.tsv
    ├── genus_sharing_across_markers.tsv
    ├── genus_sharing_barplot.pdf
    └── combined_multimarker_pcoa.pdf + .png
```

---

## Primer reference

| Marker | Primer pair | Target region | Amplicon size |
|--------|-------------|---------------|---------------|
| 16S | 515F (Parada) / 806RB (Apprill) | V4 | ~253 bp |
| 18S | TAReuk454FWD1 / TAReukREV3 | V4 | ~380 bp |
| ITS | ITS1F / ITS2 | ITS1 | ~200–500 bp |
| CO1 | mlCOIintF / jgHCO2198 (Leray) | mtCO1 | ~313 bp |
| 12S | MiFish-U-F / MiFish-U-R | 12S | ~163–185 bp |

---

## Metadata format

Optional TSV with `sample_id` as first column (matching samplesheet `sample` column):

```
sample_id    habitat    season    pH    conductivity_uS
POND_A_S1    pond       summer    7.2   312
POND_B_S2    pond       summer    7.4   298
```

Categorical columns are used for PERMANOVA grouping and plot colouring.
Numeric columns are used for RDA environmental fitting.

---

## Dependencies

All tools run via containers — no manual installation needed with `-profile docker` or `-profile singularity`.

| Tool | Version | Use |
|------|---------|-----|
| Nextflow | ≥23.04 | Workflow engine |
| FastQC | 0.12.1 | Read QC |
| MultiQC | 1.21 | QC aggregation |
| Cutadapt | 4.6 | Primer trimming |
| DADA2 | 1.30 | Denoising, taxonomy |
| VSEARCH | 2.27 | Chimera detection |
| phyloseq | ≥1.44 | Ecological analysis |
| vegan | ≥2.6 | Community ecology |
| R | 4.3 | Statistics and plots |

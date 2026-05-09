# nf-edna-pipeline

A Nextflow DSL2 pipeline for **multi-marker eDNA metabarcoding** — from raw paired-end reads to full ecological analysis. Supports **16S, 18S, ITS, CO1, 12S, and rbcL** primer sets in a single run.

---

## Contents

- [Overview](#overview)
- [Pipeline summary](#pipeline-summary)
- [Requirements](#requirements)
- [Installation](#installation)
- [Database setup](#database-setup)
- [Quick start](#quick-start)
- [Samplesheet generation](#samplesheet-generation)
- [Samplesheet format](#samplesheet-format)
- [Metadata format](#metadata-format)
- [Parameters](#parameters)
- [Profiles](#profiles)
- [Full ecological analysis](#full-ecological-analysis)
- [Output structure](#output-structure)
- [Primer reference](#primer-reference)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Overview

This pipeline processes raw Illumina amplicon reads through quality control, primer trimming, DADA2-based denoising, chimera removal, and taxonomy assignment — then feeds the resulting ASV tables into a comprehensive ecological analysis suite.

Multiple markers can be processed in a single run. The same sample can appear in the samplesheet multiple times with different marker assignments (e.g. the same DNA extract amplified with 16S and rbcL primers).

---

## Pipeline summary

```
Raw reads (FASTQ)
│
├── FastQC               — per-sample raw read QC
├── Cutadapt             — marker-specific primer trimming
├── MultiQC              — aggregated QC report
│
├── DADA2 (per run)      — learn error models
├── DADA2 (per sample)   — denoise, merge paired reads, length filter
├── VSEARCH              — de novo chimera detection (uchime3)
├── Merge ASV tables     — cross-sample merge + final bimera removal
│
├── Taxonomy assignment  (marker-specific databases)
│   ├── 16S   →  SILVA 138.1            (DADA2 naive Bayesian)
│   ├── 18S   →  PR2 v5.0              (DADA2 naive Bayesian)
│   ├── ITS   →  UNITE v10             (DADA2 naive Bayesian)
│   ├── CO1   →  MIDORI2 (GB262)       (DADA2 naive Bayesian)
│   ├── 12S   →  MIDORI2 (GB262)       (DADA2 naive Bayesian)
│   └── rbcL  →  rbcLClassifier v1     (RDP Classifier, NCBI-trained)
│
├── Basic ecology        — diversity, PCoA, NMDS, barplots
│
└── Full ecology suite   — 9 modules (optional, runs after basic)
    ├── Alpha diversity      rarefaction, richness, evenness, statistics
    ├── Beta diversity       6 distance metrics, PERMANOVA, PERMDISP, ANOSIM, MRPP
    ├── Ordination           PCoA, NMDS, PCA-CLR, RDA, db-RDA, CCA
    ├── Differential abund.  DESeq2 + ALDEx2 (all pairwise)
    ├── Co-occurrence net.   Spearman network, modularity, hub taxa
    ├── Indicator species    IndVal, SIMPER, core microbiome
    ├── Environmental drivers envfit, Mantel, variance partitioning, Procrustes
    ├── Cross-marker         inter-marker Procrustes, genus sharing, combined PCoA
    └── HTML report          per-marker integrated R Markdown report
```

---

## Requirements

| Software | Minimum version | Notes |
|----------|----------------|-------|
| [Nextflow](https://nextflow.io) | 23.04.0 | `curl -s get.nextflow.io \| bash` |
| [Docker](https://docs.docker.com/get-docker/) **or** [Singularity](https://sylabs.io/guides/3.0/user-guide/) | Docker ≥ 20 / Singularity ≥ 3.7 | One or the other required for containers |
| Java | 11–21 | Required by Nextflow |
| Python | 3.8+ | Required for `bin/make_samplesheet.py` (standard library only) |

All bioinformatics tools (FastQC, Cutadapt, DADA2, VSEARCH, RDP Classifier, etc.) are pulled automatically as containers — no manual tool installation is needed.

---

## Installation

**1. Install Nextflow**

```bash
# Install to ~/bin (or any directory on your PATH)
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/
```

Verify:
```bash
nextflow -version
```

**2. Clone this repository**

```bash
git clone git@github.com:UKCEH-MolecularEcology/nf-edna-pipeline.git
cd nf-edna-pipeline
```

**3. Install a container engine**

*Docker (local workstations):*
```bash
# Ubuntu/Debian
sudo apt-get install -y docker.io
sudo usermod -aG docker $USER   # log out and back in after this
```

*Singularity (HPC clusters — usually pre-installed):*
```bash
singularity --version   # check it's available
```

---

## Database setup

Taxonomy assignment requires reference databases. A helper script downloads them automatically:

```bash
bash assets/download_databases.sh databases/
```

This downloads:
- **SILVA 138.1** (16S) from Zenodo
- **PR2 v5.0** (18S) from GitHub
- **MIDORI2 GB262** (CO1, 12S) from the MIDORI website
- **rbcLClassifier v1** (rbcL) from GitHub — RDP Classifier trained on NCBI plant rbcL sequences (Kress & Porter 2020); citable archive at [doi:10.5281/zenodo.4741459](https://doi.org/10.5281/zenodo.4741459)

> **UNITE (ITS)** requires manual download due to licensing. Go to [unite.ut.ee/repository.php](https://unite.ut.ee/repository.php), download the *General FASTA release — developer s_all version*, and place it in `databases/`. Update the path in `nextflow.config` accordingly.

After downloading, verify the paths in `nextflow.config` under `params.databases` match your local files:

```groovy
databases {
    '16S'  { path = "${projectDir}/databases/silva_nr99_v138.1_train_set.fa.gz" }
    '18S'  { path = "${projectDir}/databases/pr2_version_5.0.0_SSU_dada2.fasta.gz" }
    'ITS'  { path = "${projectDir}/databases/sh_general_release_dynamic_s_all_19.02.2024.fasta.gz" }
    'CO1'  { path = "${projectDir}/databases/MIDORI2_LONGEST_NUC_GB262_CO1_DADA2.fasta.gz" }
    '12S'  { path = "${projectDir}/databases/MIDORI2_LONGEST_NUC_GB262_12S_DADA2.fasta.gz" }
    'RBCL' { path = "${projectDir}/databases/rbcLv1_trained/mydata_trained" }
}
```

---

## Quick start

There are two ways to provide input reads. Use `--input` (samplesheet) or `--fastq_dir` (directory) — not both.

**Option A — auto-detect from a FASTQ directory (no samplesheet needed):**

```bash
# Markers are inferred from filenames; a samplesheet is written to results/
nextflow run main.nf \
    --fastq_dir /path/to/fastqs/ \
    --cores 16 \
    --outdir results/ \
    -profile docker

# With metadata for ecology analysis on an HPC cluster
nextflow run main.nf \
    --fastq_dir /path/to/fastqs/ \
    --metadata metadata.tsv \
    --ecology_group_var habitat \
    --cores 64 \
    --outdir results/ \
    -profile slurm,singularity \
    -resume
```

**Option B — explicit samplesheet:**

```bash
# Run 16S with Docker
nextflow run main.nf \
    --input samplesheet.csv \
    --markers 16S \
    --cores 8 \
    --outdir results/ \
    -profile docker

# Run all six markers on a SLURM cluster
nextflow run main.nf \
    --input samplesheet.csv \
    --markers 16S,18S,ITS,CO1,12S,RBCL \
    --metadata metadata.tsv \
    --ecology_group_var habitat \
    --cores 64 \
    --outdir results/ \
    -profile slurm,singularity \
    -resume
```

Use `-resume` to restart from the last successful step after a failure or parameter change.

> **Without metadata:** All quality, denoising, and taxonomy steps run normally. Ecology analyses run in unsupervised mode — alpha diversity, ordination, barplots, and co-occurrence networks are produced, but group-dependent analyses (PERMANOVA, differential abundance, IndVal, envfit) are automatically skipped.

---

## Samplesheet generation

The pipeline can consume a FASTQ directory directly via `--fastq_dir` (see [Quick start](#quick-start)), which auto-detects markers and writes the samplesheet to `results/samplesheet_detected.csv` for reference.

Alternatively, use the helper script to generate a samplesheet ahead of time — useful for reviewing or editing before running the pipeline.

**Expected filename format** (standard Illumina output):
```
{sample_id}-{marker}_{S-index}_{R1|R2}_001.fastq.gz
```

**Examples:**
```
H4D-EA100_2024-18S_S244_R1_001.fastq.gz
K9G-EA71_2024-rbcL_S503_R1_001.fastq.gz
P8B-Blank2-SW-rbcL_S762_R2_001.fastq.gz
```

The script detects the marker from the last `-`-delimited token before the Illumina suffix and recognises: `16S`, `18S`, `ITS`, `CO1`/`COI`, `12S`, `RBCL`.

**Generate samplesheet only:**
```bash
python bin/make_samplesheet.py /path/to/fastqs/ -o samplesheet.csv
```

**Generate samplesheet and filter metadata to matched samples:**
```bash
python bin/make_samplesheet.py /path/to/fastqs/ \
    -o samplesheet.csv \
    --metadata metadata.tsv
```

This produces:
- `samplesheet.csv` — ready to pass to `--input`
- `samplesheet_metadata.tsv` — filtered metadata containing only samples found in the FASTQ directory, ready to pass to `--metadata`
- A warning for any samples present in one file but missing from the other
- A ready-to-run `nextflow run` command printed to stdout

**Options:**
```
positional:  fastq_dir          Directory containing *.fastq.gz files
-o/--output  samplesheet.csv    Output CSV path (default: samplesheet.csv)
-m/--metadata                   Metadata TSV/CSV to merge
--metadata-out                  Custom path for filtered metadata output
--absolute                      Write absolute FASTQ paths
```

---

## Samplesheet format

Create a comma-separated file with these columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample identifier — no spaces or special characters |
| `fastq_1` | Yes | Absolute path to R1 FASTQ file (gzipped) |
| `fastq_2` | No | Absolute path to R2 FASTQ file (gzipped). Omit for single-end |
| `marker` | Yes | One of: `16S`, `18S`, `ITS`, `CO1`, `12S`, `RBCL` |
| `run` | No | Sequencing run ID. Used to group error-model learning. Defaults to `run1` |

**Example (`samplesheet.csv`):**

```csv
sample,fastq_1,fastq_2,marker,run
POND_A,/data/POND_A_R1.fastq.gz,/data/POND_A_R2.fastq.gz,16S,run1
POND_A,/data/POND_A_R1.fastq.gz,/data/POND_A_R2.fastq.gz,RBCL,run1
RIVER_B,/data/RIVER_B_R1.fastq.gz,/data/RIVER_B_R2.fastq.gz,16S,run1
RIVER_B,/data/RIVER_B_R1.fastq.gz,/data/RIVER_B_R2.fastq.gz,CO1,run1
MEADOW_C,/data/MEADOW_C_R1.fastq.gz,/data/MEADOW_C_R2.fastq.gz,RBCL,run2
NEG_CTRL,/data/NEG_R1.fastq.gz,/data/NEG_R2.fastq.gz,16S,run1
```

> **Note:** The same physical FASTQ files can appear multiple times (one row per marker per sample). The pipeline treats each row as an independent amplification of the same extract.

> **Sequencing runs:** If your samples span multiple sequencing runs, use the `run` column. DADA2 error models are learned per run — mixing runs without the `run` column will reduce accuracy.

---

## Metadata format

An optional tab-separated metadata file enriches the ecological analysis with groupings, environmental variables, and statistical comparisons.

The first column must be named `sample_id` and match the `sample` values in the samplesheet. The first categorical column is used as the primary grouping variable by default (override with `--ecology_group_var`).

```tsv
sample_id	habitat	season	pH	conductivity_uS	dissolved_O2_mg_L
POND_A	pond	summer	7.2	312	8.1
RIVER_B	river	summer	7.8	485	9.2
MEADOW_C	meadow	summer	6.5	NA	NA
```

- **Categorical columns** — used for PERMANOVA grouping, ANOSIM, IndVal, boxplot colouring, differential abundance
- **Numeric columns** — used for RDA/CCA environmental fitting, Mantel tests, variance partitioning, envfit vectors, alpha diversity correlations

Use `bin/make_samplesheet.py --metadata` to automatically filter your metadata file to only the samples present in your FASTQ directory before passing it to the pipeline.

---

## Parameters

### Core parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | Samplesheet CSV. Mutually exclusive with `--fastq_dir` |
| `--fastq_dir` | — | Directory of FASTQ files. Markers auto-detected from filenames. Mutually exclusive with `--input` |
| `--outdir` | `results` | Output directory |
| `--markers` | `16S` | Comma-separated markers to process. Ignored when using `--fastq_dir` |
| `--metadata` | null | Path to sample metadata TSV. If omitted, ecology runs in unsupervised mode |
| `--cores` | `4` | Number of CPU threads available to the pipeline. All process labels scale from this value (see table below) |
| `--single_end` | `false` | Set to `true` for single-end libraries |

**CPU scaling by `--cores`:**

| `--cores` | `process_low` (Cutadapt, VSEARCH) | `process_medium` (DADA2 error models) | `process_high` (DADA2 denoising, taxonomy, ecology) |
|-----------|----------------------------------|--------------------------------------|-----------------------------------------------------|
| 4 | 2 | 4 | 4 |
| 8 | 2 | 4 | 8 |
| 16 | 2 | 4 | 16 |
| 32 | 4 | 8 | 32 |
| 64 | 8 | 16 | 64 |

### DADA2 parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dada2_pool` | `false` | Pooling strategy: `false` (independent), `true` (full pooling), or `pseudo` (recommended for rare species detection) |
| `--dada2_chimera` | `consensus` | Chimera removal method for merged table: `consensus`, `pooled`, or `per-sample` |

### RDP Classifier parameters (rbcL)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--rdp_mem` | `8g` | JVM heap size for RDP Classifier (e.g. `8g`, `16g`) |
| `--rdp_bootstrap` | `0.70` | Minimum bootstrap confidence (0–1) to retain a rank assignment. 0.70 gives ~90% accuracy at genus level for ≥400 bp reads |

### Primer parameters

Default primer sequences are defined in `nextflow.config` under `params.primers` for all six markers. Override any primer or quality parameter at the command line:

```bash
nextflow run main.nf \
    --markers RBCL \
    --primers.RBCL.fwd ATGTCACCACAAACAGAGACTAAAGC \
    --primers.RBCL.rev GTAAAATCAAGTCCACCRCG \
    --primers.RBCL.min_length 400 \
    ...
```

Or edit `nextflow.config` directly for permanent changes.

### Ecology parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_ecology` | `true` | Run basic ecological analysis (diversity, PCoA, barplots) |
| `--run_full_ecology` | `true` | Run the full 9-module ecological analysis suite |
| `--ecology_group_var` | null | Metadata column to use as primary group variable. Auto-detects first categorical column if unset |
| `--ecology_min_prevalence` | `0.3` | Network analysis: minimum fraction of samples an ASV must appear in |
| `--ecology_cor_cutoff` | `0.6` | Network analysis: minimum Spearman \|r\| threshold for a co-occurrence edge |

### Resource limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `32` | Hard upper limit on CPUs per process |
| `--max_memory` | `128.GB` | Hard upper limit on memory per process |
| `--max_time` | `240.h` | Hard upper limit on walltime per process |

---

## Profiles

Select a profile with `-profile <name>`. Combine multiple profiles with commas (e.g. `-profile slurm,singularity`).

| Profile | Description |
|---------|-------------|
| `standard` | Local execution, no containers. Requires all tools installed manually |
| `docker` | Local execution with Docker containers |
| `singularity` | Container execution via Singularity (recommended for HPC) |
| `slurm` | SLURM HPC scheduler + Singularity. Adjust queue names in `nextflow.config` if needed |
| `pbs` | PBS/Torque HPC scheduler + Singularity |
| `test` | Minimal test run with bundled test data and reduced resources |

**Example — SLURM cluster:**

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --markers 16S,RBCL \
    --metadata metadata.tsv \
    --outdir results/ \
    -profile slurm,singularity \
    -resume
```

---

## Full ecological analysis

The full ecology suite runs automatically after the main pipeline when `--run_full_ecology true` (default). It can also be run **standalone** against results from a previous pipeline run — useful for re-running analysis with different parameters without re-processing reads.

### Standalone usage

```bash
# Point at a previous results directory
nextflow run ecology_pipeline.nf \
    --results_dir results/ \
    --markers 16S,RBCL \
    --metadata metadata.tsv \
    --ecology_group_var habitat \
    --outdir full_ecology/ \
    -profile singularity

# Or supply files directly (single marker)
nextflow run ecology_pipeline.nf \
    --asv_table results/asv_tables/RBCL/RBCL.merged_asv_table.tsv \
    --taxonomy  results/taxonomy/RBCL/SAMPLE_RBCL.taxonomy.tsv \
    --marker    RBCL \
    --metadata  metadata.tsv \
    --outdir    full_ecology/ \
    -profile docker
```

### What each module produces

| Module | Key outputs |
|--------|------------|
| **Alpha diversity** | Richness/evenness/diversity metrics TSV, boxplots, Kruskal-Wallis + Dunn's post-hoc, iNEXT rarefaction curves, rank-abundance curves, occupancy-abundance plot |
| **Beta diversity** | 6 distance matrix TSVs, distance heatmap, PERMANOVA (single + multi-factor), pairwise PERMANOVA (BH-adjusted), PERMDISP, ANOSIM, MRPP, multivariate test summary |
| **Ordination** | PCoA (Bray-Curtis + Aitchison), NMDS + Shepard diagram, PCA-CLR biplot with taxon loadings, RDA with forward variable selection + axis ANOVA, db-RDA, CCA, 4-panel summary PDF |
| **Differential abundance** | DESeq2 all pairwise TSVs + volcano plots + MA plot + significant-taxa heatmap; ALDEx2 results + effect plot + volcano |
| **Co-occurrence network** | Spearman correlation matrix, edge list, network statistics, node metrics, hub/keystone taxa, Louvain modularity, ggraph network plot |
| **Indicator species** | IndVal significant indicators + barplot, SIMPER tables for all group pairs + contribution plot, core microbiome matrix (9 threshold combinations) |
| **Environmental drivers** | envfit vectors on NMDS, Mantel test per variable, variance partitioning diagram, Procrustes (community vs env PCA) plot, alpha diversity × env Spearman heatmap |
| **Cross-marker** | Procrustes correlation between all marker pairs, genus sharing matrix, combined multi-marker PCoA |
| **HTML report** | Per-marker R Markdown report embedding all tables and figures |

---

## Output structure

```
results/
├── pipeline_info/                     Execution report, timeline, trace, DAG
├── fastqc/{MARKER}/                   Per-sample FastQC HTML + ZIP
├── multiqc/                           Aggregated QC HTML report
├── trimmed/{MARKER}/                  Cutadapt logs
├── dada2/{MARKER}/
│   ├── error_models/                  Error model RDS + diagnostic plots
│   └── asv_tables/                    Per-sample ASV table RDS, FASTA, read stats
├── chimera_check/{MARKER}/            VSEARCH chimera output
├── asv_tables/{MARKER}/
│   ├── {MARKER}.merged_asv_table.tsv  ← Final ASV count table (ASVs × samples)
│   ├── {MARKER}.merged_asv_table.rds  Phyloseq-ready RDS
│   └── {MARKER}.read_tracking.tsv     Reads surviving each step
├── taxonomy/{MARKER}/
│   └── {SAMPLE}_{MARKER}.taxonomy.tsv Taxonomy assignments with bootstrap scores
├── ecology/{MARKER}/                  Basic ecology (diversity, ordination, barplots)
└── full_ecology/
    ├── {MARKER}/
    │   ├── 00_report/                 {MARKER}_ecological_report.html  ← Start here
    │   ├── 01_alpha_diversity/
    │   ├── 02_beta_diversity/
    │   ├── 03_ordination/
    │   ├── 04_differential_abundance/
    │   ├── 05_co_occurrence_network/
    │   ├── 06_indicator_species/
    │   └── 07_env_drivers/
    └── cross_marker/                  Inter-marker comparisons
```

> **Tip:** Start your exploration with the HTML report at `full_ecology/{MARKER}/00_report/{MARKER}_ecological_report.html`. It links to all tables and figures.

---

## Primer reference

Default primers used by the pipeline. To use different primers, override `params.primers` in `nextflow.config` or via `--primers.{MARKER}.fwd` / `--primers.{MARKER}.rev` on the command line.

| Marker | Forward primer | Reverse primer | Target region | Expected amplicon |
|--------|---------------|----------------|---------------|-------------------|
| 16S | 515F (Parada): `GTGYCAGCMGCCGCGGTAA` | 806RB (Apprill): `GGACTACNVGGGTWTCTAAT` | V4 | ~253 bp |
| 18S | TAReuk454FWD1: `CCAGCASCYGCGGTAATTCC` | TAReukREV3: `ACTTTCGTTCTTGATYRA` | V4 | ~380 bp |
| ITS | ITS1F: `CTTGGTCATTTAGAGGAAGTAA` | ITS2: `GCTGCGTTCTTCATCGATGC` | ITS1 | 200–500 bp |
| CO1 | mlCOIintF (Leray): `GGWACWGGWTGAACWGTWTAYCCYCC` | jgHCO2198: `TAIACYTCIGGRTGICCRAARAAYCA` | mtCO1 | ~313 bp |
| 12S | MiFish-U-F: `GTCGGTAAAACTCGTGCCAGC` | MiFish-U-R: `CATAGTGGGGTATCTAATCCCAGTTTG` | 12S rRNA | 163–185 bp |
| rbcL | rbcLa-F (Kress & Erickson 2007): `ATGTCACCACAAACAGAGACTAAAGC` | rbcLa-R: `GTAAAATCAAGTCCACCRCG` | rbcLa | ~550 bp |

---

## Troubleshooting

**Pipeline fails at DADA2 with "not enough sequences"**
: Increase `max_ee_f` / `max_ee_r` (e.g. to 3–5) or reduce `trunc_len_f` / `trunc_len_r`. Check the Cutadapt log — low trimmed read counts indicate primers were not found.

**DADA2 error model produces warnings about insufficient reads**
: Each sequencing run needs ≥ ~100,000 reads across all samples to learn a reliable error model. If you have very few samples, try `--dada2_pool pseudo`.

**rbcL classification is slow or runs out of memory**
: Increase `--rdp_mem` (e.g. `--rdp_mem 16g`). The RDP Classifier loads the entire trained model into JVM heap; 8 GB is usually sufficient for the rbcLClassifier v1.

**rbcL taxonomy assignments are mostly blank**
: Lower `--rdp_bootstrap` (e.g. `0.50`). The default 0.70 is conservative; lower values retain more assignments at the cost of reduced accuracy at finer ranks. Assignments are only possible for reads ≥ 400 bp after trimming.

**Taxonomy database not found**
: Check the path in `nextflow.config` matches the exact filename/directory downloaded. Paths support `${projectDir}` as a shorthand for the pipeline directory. For rbcL, the path should point to the `mydata_trained/` directory extracted from `rbcLv1_trained.tar.gz`.

**OutOfMemoryError in DADA2 or taxonomy**
: Increase `--max_memory` (e.g. `--max_memory 256.GB`) or set `process.memory` directly in a custom config.

**SLURM jobs queuing but not running**
: Check the queue name in `nextflow.config` under `profiles.slurm`. Adjust `process.queue` to match your cluster's partition names.

**Resuming a pipeline run**
: Always use `-resume` when re-running after a failure. Nextflow caches completed process outputs in `work/` and skips them.

**Work directory getting large**
: Clean up completed runs with `nextflow clean -f` or delete `work/` manually. Only do this once you no longer need to `-resume`.

---

## Citation

If you use this pipeline in your research, please cite:

- **DADA2:** Callahan BJ *et al.* (2016) DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581–583. https://doi.org/10.1038/nmeth.3869
- **VSEARCH:** Rognes T *et al.* (2016) VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584. https://doi.org/10.7717/peerj.2584
- **SILVA:** Quast C *et al.* (2013) The SILVA ribosomal RNA gene database project. *Nucleic Acids Research* 41:D590–D596. https://doi.org/10.1093/nar/gks1219
- **PR2:** Guillou L *et al.* (2013) The Protist Ribosomal Reference database (PR2). *Nucleic Acids Research* 41:D597–D604. https://doi.org/10.1093/nar/gks1160
- **UNITE:** UNITE Community (2023) UNITE general FASTA release. https://doi.org/10.15156/BIO/2938065
- **MIDORI2:** Machida RJ *et al.* (2017) MIDORI series of curated sequence reference databases. *Scientific Data* 4:170027. https://doi.org/10.1038/sdata.2017.27
- **rbcLClassifier:** Porter TM & Hajibabaei M (2022) rbcLClassifier: a trained RDP Classifier for rbcL (v1.0). Zenodo. https://doi.org/10.5281/zenodo.4741459
- **RDP Classifier:** Wang Q *et al.* (2007) Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and Environmental Microbiology* 73:5261–5267. https://doi.org/10.1128/AEM.00062-07
- **rbcL primers:** Kress WJ & Erickson DL (2007) A two-locus global DNA barcode for land plants. *PLOS ONE* 2:e508. https://doi.org/10.1371/journal.pone.0000508
- **phyloseq:** McMurdie PJ & Holmes S (2013) phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE* 8:e61217. https://doi.org/10.1371/journal.pone.0061217
- **vegan:** Oksanen J *et al.* (2022) vegan: Community Ecology Package. R package. https://CRAN.R-project.org/package=vegan
- **DESeq2:** Love MI *et al.* (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550. https://doi.org/10.1186/s13059-014-0550-8
- **ALDEx2:** Fernandes AD *et al.* (2014) Unifying the analysis of high-throughput sequencing datasets. *PLOS Computational Biology* 10:e1003531. https://doi.org/10.1371/journal.pcbi.1003531
- **iNEXT:** Hsieh TC *et al.* (2016) iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers). *Methods in Ecology and Evolution* 7:1451–1456. https://doi.org/10.1111/2041-210X.12613
- **indicspecies:** De Cáceres M & Legendre P (2009) Associations between species and groups of sites. *Ecology* 90:3566–3574. https://doi.org/10.1890/08-1823.1
- **Nextflow:** Di Tommaso P *et al.* (2017) Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35:316–319. https://doi.org/10.1038/nbt.3820

---

*Developed by the Molecular Ecology group at the UK Centre for Ecology & Hydrology (UKCEH).*

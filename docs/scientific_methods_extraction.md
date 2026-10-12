# Scientific Methods Extraction Report: nf-core/methylarray

> Auto-generated documentation for writing a peer-reviewed scientific methods paper.
> Extracted from pipeline source code, configuration files, and R scripts.
> Last updated: 2026-03-26

---

## SECTION 1 — PIPELINE IDENTITY

**Pipeline name:** nf-core/methylarray

**Version:** 1.0.0dev (CHANGELOG entry: v1.0.0dev, 2026-03-15)

**DOI:** [NOT FOUND — placeholder `10.5281/zenodo.XXXXXXX` in README; no real DOI assigned yet]

**nf-core URL:** https://nf-co.re/methylarray

**GitHub URL:** https://github.com/nf-core/methylarray

**Slack channel:** https://nfcore.slack.com/channels/methylarray

**Authors and affiliations** (from `nextflow.config` manifest):

| Name            | Affiliation           | Email                         | GitHub                           | Role               | ORCID   |
| --------------- | --------------------- | ----------------------------- | -------------------------------- | ------------------ | ------- |
| Adam Schumacher | Karolinska Institutet | adam.schumacher.pro@gmail.com | https://github.com/schumz        | author, maintainer | [blank] |
| Ghada Nouairia  | Karolinska Institutet | ghada.nouairia@ki.se          | https://github.com/GhadaNOUAIRIA | author             | [blank] |

**License:** MIT (Copyright: The nf-core/methylarray team)

**One-line description:** "Pipeline for analysing Illumina EPIC v2.0 DNA methylation IDAT data. It includes sample and probe QC (detection p-values, signal intensities, control probes) with filtering of low-quality samples/probes, followed by background correction and normalization, and downstream differential methylation analysis to identify DMPs and DMRs."

**Date of last release:** 2026-03-15 (v1.0.0dev — initial development release)

**nf-core template version:** 3.5.2

**Required Nextflow version:** ≥ 25.04.0 (hard requirement: `!>=25.04.0`)

---

## SECTION 2 — BIOLOGICAL / SCIENTIFIC PURPOSE

**Biological question addressed:**
The pipeline addresses epigenome-wide association study (EWAS) of DNA methylation in human samples — identifying differentially methylated CpG positions (DMPs) and regions (DMRs) between disease groups and controls. It is designed specifically for blood-based epigenetic biomarker discovery where cell-type heterogeneity is a major confounder.

**Type of input data:**

- Raw IDAT files from Illumina EPICv2 (~900K probes) DNA methylation arrays (both green and red channel raw fluorescence intensity files)
- A tab-separated samplesheet linking sample names to Sentrix array coordinates (Sample_Name, Sentrix_ID, Sentrix_Position)
- An optional metadata CSV with clinical/phenotype information per sample

**Expected outputs:**

- QC reports (raw intensity plots, control probe strip plots, SNP identity heatmaps, sex concordance plots, beta-value density plots)
- Normalized beta-value and M-value matrices (RDS format, optionally gzip CSV)
- Blood cell-type composition estimates (CD8T, CD4T, NK, Bcell, Mono, Neu, Gran)
- Covariate-PC association tables and ChAMP SVD outputs
- DMP tables per pairwise group contrast (limma, three batch-correction strategies)
- DMR tables per pairwise group contrast (DMRcate, three batch-correction strategies)
- Optional consensus DMRs from both DMRcate and coMethDMR
- Comprehensive pipeline execution reports (timeline, trace, DAG, software versions)

**Organism(s):** Human (Homo sapiens), GRCh38/hg38 (default); hg19 supported as legacy option.

**Array technology:** Illumina Infinium MethylationEPIC v2.0 (EPICv2, ~900K probes), with backward-compatible support for EPICv1 (850K). Default annotation: `IlluminaHumanMethylationEPICv2` (`20a1.hg38`).

**Experimental context:** Designed for blood-based EWAS in cohort studies with a case-control design (multiple disease groups vs. CONTROL). Supports multi-group comparisons with all pairwise contrasts computed automatically.

---

## SECTION 3 — WORKFLOW ARCHITECTURE

### Workflow files and their purpose

| File                                                           | Purpose                                                                                                                                                  |
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `main.nf`                                                      | Entry point: initialises the pipeline via `PIPELINE_INITIALISATION` subworkflow, calls `NFCORE_METHYLARRAY` wrapper workflow, then `PIPELINE_COMPLETION` |
| `workflows/methylarray.nf`                                     | Main analysis workflow `METHYLARRAY`: imports all 16 local modules and orchestrates the full pipeline from IDAT import to DMR calling                    |
| `subworkflows/local/utils_nfcore_methylarray_pipeline/main.nf` | `PIPELINE_INITIALISATION` (param validation, channel creation) and `PIPELINE_COMPLETION` (email/hook notifications)                                      |
| `subworkflows/nf-core/utils_nextflow_pipeline`                 | nf-core utility for Nextflow-level operations (version display, parameter dump)                                                                          |
| `subworkflows/nf-core/utils_nfcore_pipeline`                   | nf-core utility for nf-core-specific operations (config checks, completion email)                                                                        |

### Step-by-step execution flow

```
INPUT
  ↓
1.  PIPELINE_INITIALISATION   — Param validation, channel creation
2.  LOG_VERSIONS              — R/Bioconductor package version logging
3.  IMPORT_IDAT               — IDAT loading → RGChannelSet (rg0.rds), targets.tsv
4.  RAW_INTENSITIES_QC        — Median intensity QC (flag low-quality samples)
5.  SESAME_POOBAH_QC          — SeSAMe QCDPB normalization + pOOBAH masking + sample/probe removal → beta_norm.rds
6.  QC_REPORT_CONTROLS        — minfi control probe PDFs and strip plots
7.  SNP_HEATMAP               — 59 SNP probe correlation heatmap (identity check)
8.  DENSITY_PLOTS             — Raw vs SeSAMe beta density overlay plots
9.  SEX_QC                    — Predicted sex vs reported sex concordance
10. FILTER_XY_NONCG_BEADS     — Bead count / sex chromosome / non-CpG probe removal → beta3.rds
11. CELL_COMP*                — Houseman blood deconvolution (IDOL CpGs)
12. SPLIT_COLLAPSE            — SNP/CH removal; EPICv2 replicate collapse (DMP branch) → beta_clean.rds
                              — SNP removal; no collapse (DMR branch) → beta_for_DMR_uncollapsed.rds
13. FINAL_EXPORTS             — M-value computation; probe/sample manifests
14. ALIGN_META                — Join beta/M matrices with metadata, predicted sex, cell composition
15. COVARIATE_PCA*            — PCA on top 50k variable CpGs; PC-covariate association; ChAMP SVD*
16. DMP_LIMMA*                — limma DMP analysis (3 modes: norm_only, plate_in_model, combat+SVA)
17. DMR_DMRCATE*              — DMRcate DMR analysis (3 modes); optional coMethDMR + consensus*
OUTPUT
```

_\* = conditional step_

### Conditional steps

| Step                  | Condition                                                                   |
| --------------------- | --------------------------------------------------------------------------- |
| CELL_COMP             | `params.do_estimate_cellcomp == true` (Groovy `if` in workflow)             |
| COVARIATE_PCA         | `ext.when = { params.run_covariate_pca }` in modules.config                 |
| DMP_LIMMA             | `ext.when = { params.run_dmp }` in modules.config                           |
| DMR_DMRCATE           | `ext.when = { params.run_dmr }` in modules.config                           |
| coMethDMR + consensus | `params.run_cometh == true` (checked inside dmr_dmrcate.R)                  |
| ChAMP SVD             | `params.run_champ_svd == true` (passed as `--run-champ` to covariate_pca.R) |
| `rg0.rds` saved       | `params.save_intermediates == true`                                         |
| `beta_norm.rds` saved | `params.save_intermediates == true`                                         |
| `beta3.rds` saved     | `params.save_intermediates == true`                                         |
| CSV matrix exports    | `params.export_csv_matrices == true`                                        |

---

## SECTION 4 — MODULES AND TOOLS

| Tool                                             | Version                  | Container / Conda source                                                 | Pipeline Step                                                                                                             | Role                                                                                                                                                                                                             | Key Parameters                                                                                                |
| ------------------------------------------------ | ------------------------ | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| **minfi**                                        | 1.48.0                   | `quay.io/nf-core/methylarray:1.0.0dev` / `bioconductor-minfi=1.48.0`     | IMPORT_IDAT, RAW_INTENSITIES_QC, QC_REPORT_CONTROLS, SNP_HEATMAP, DENSITY_PLOTS, SEX_QC, FILTER_XY_NONCG_BEADS, CELL_COMP | IDAT loading (`read.metharray.exp`), raw QC (`preprocessRaw`, `getQC`), control probe plots, SNP beta extraction, density plots, sex prediction (`getSex`, `preprocessNoob`), bead count filtering (`getNBeads`) | `extended=TRUE`; annotation forced to EPICv2                                                                  |
| **SeSAMe**                                       | 1.18.0                   | same / `bioconductor-sesame=1.18.0`                                      | SESAME_POOBAH_QC, SPLIT_COLLAPSE                                                                                          | Preprocessing pipeline (`openSesame`), pOOBAH masking (`getBetas`), EPICv2 replicate collapsing (`betasCollapseToPfx`)                                                                                           | `--sesame-prep QCDPB`, `--p-det 0.05`, `BiocParallel::SerialParam()`                                          |
| **sesameData**                                   | 1.18.0                   | same / `bioconductor-sesamedata=1.18.0`                                  | SESAME_POOBAH_QC                                                                                                          | Data cache for SeSAMe quality masks and reference data (`sesameDataCache`)                                                                                                                                       | Local hub mode enforced                                                                                       |
| **IlluminaHumanMethylationEPICv2manifest**       | 1.0.0                    | same / `bioconductor-illuminahumanmethylationepicv2manifest=1.0.0`       | IMPORT_IDAT, FILTER_XY_NONCG_BEADS                                                                                        | EPICv2 array manifest for probe definitions and channel assignment                                                                                                                                               | annotation=`20a1.hg38`                                                                                        |
| **IlluminaHumanMethylationEPICv2anno.20a1.hg38** | 1.0.0                    | same / `bioconductor-illuminahumanmethylationepicv2anno.20a1.hg38=1.0.0` | IMPORT_IDAT, FILTER_XY_NONCG_BEADS                                                                                        | Genomic annotation for EPICv2 probes (chromosome, position) — used for XY filtering                                                                                                                              | `--annotation 20a1.hg38`                                                                                      |
| **limma**                                        | 3.58.0                   | same / `bioconductor-limma=3.58.0`                                       | DMP_LIMMA, DMR_DMRCATE                                                                                                    | Linear model fitting (`lmFit`), empirical Bayes moderation (`eBayes`), contrast definition, results extraction (`topTable`)                                                                                      | `trend=TRUE`, `robust=TRUE`; BH p-value adjustment                                                            |
| **DMRcate**                                      | 2.16.0                   | same / `bioconductor-dmrcate=2.16.0`                                     | SPLIT_COLLAPSE, DMR_DMRCATE                                                                                               | SNP/cross-hybridisation removal (`rmSNPandCH`), CpG annotation (`cpg.annotate`), DMR calling (`dmrcate`), range extraction (`extractRanges`)                                                                     | `lambda=1000`, `C=2`, `epicv2Filter=posreps_strategy`, `epicv2Remap=TRUE`, `fdr=0.05`, `arraytype="EPICv2"`   |
| **FlowSorted.Blood.EPIC**                        | 2.4.0                    | same / `bioconductor-flowsorted.blood.epic=2.4.0`                        | CELL_COMP                                                                                                                 | Reference-based blood cell deconvolution using IDOL optimized CpGs (`IDOLOptimizedCpGs`, `projectCellType_CP`)                                                                                                   | `nonnegative=TRUE`, `lessThanOne=FALSE`; min 200 IDOL CpG overlap required                                    |
| **ChAMP**                                        | 2.32.0                   | same / `bioconductor-champ=2.32.0`                                       | COVARIATE_PCA                                                                                                             | SVD-based covariate association analysis (`champ.SVD()`)                                                                                                                                                         | 1% NA rate threshold before SVD; local hub mode enforced                                                      |
| **sva**                                          | 3.50.0                   | same / `bioconductor-sva=3.50.0`                                         | DMP_LIMMA, DMR_DMRCATE                                                                                                    | ComBat batch correction (`ComBat`), surrogate variable estimation (`num.sv`, `sva`)                                                                                                                              | `par.prior=TRUE`, `prior.plots=FALSE`; "be" and "leek" k-selection methods; backoff loop for full-rank design |
| **missMethyl**                                   | 1.36.0                   | same / `bioconductor-missmethyl=1.36.0`                                  | DMP_LIMMA (installed, not yet called)                                                                                     | Gene set testing for methylation data — available for future GO/pathway enrichment                                                                                                                               | [reserved for future use]                                                                                     |
| **coMethDMR**                                    | [NOT FOUND — not pinned] | same                                                                     | DMR_DMRCATE (`--run_cometh true`)                                                                                         | Co-methylation DMR calling (`CoMethAllRegions`); consensus with DMRcate via `GenomicRanges::findOverlaps`                                                                                                        | `betaToM=FALSE`, `BiocParallel::SerialParam()`, FDR < 0.05                                                    |
| **pheatmap**                                     | 1.0.12                   | same                                                                     | SNP_HEATMAP                                                                                                               | Clustered heatmap of SNP probe correlation matrix                                                                                                                                                                | Default hierarchical clustering                                                                               |
| **ggplot2**                                      | 3.5.1                    | same / `r-ggplot2=3.5.1`                                                 | DMP_LIMMA, DMR_DMRCATE, COVARIATE_PCA                                                                                     | Plotting infrastructure                                                                                                                                                                                          | —                                                                                                             |
| **dplyr**                                        | 1.1.4                    | same / `r-dplyr=1.1.4`                                                   | ALIGN_META, COVARIATE_PCA, DMP_LIMMA, DMR_DMRCATE                                                                         | Data manipulation (filtering, joining, selecting)                                                                                                                                                                | —                                                                                                             |
| **readr**                                        | 2.1.5                    | same / `r-readr=2.1.5`                                                   | All R modules                                                                                                             | CSV/TSV reading and writing                                                                                                                                                                                      | `show_col_types=FALSE` throughout                                                                             |
| **tibble**                                       | 3.2.1                    | same / `r-tibble=3.2.1`                                                  | DMP_LIMMA, ALIGN_META                                                                                                     | Row-name handling (`rownames_to_column`, `as_tibble`)                                                                                                                                                            | —                                                                                                             |
| **matrixStats**                                  | 1.3.0                    | same / `r-matrixstats=1.3.0`                                             | COVARIATE_PCA                                                                                                             | Row variance computation (`rowVars`) for variable CpG selection                                                                                                                                                  | —                                                                                                             |
| **optparse**                                     | 1.7.5                    | same / `r-optparse=1.7.5`                                                | All R modules                                                                                                             | Command-line argument parsing                                                                                                                                                                                    | —                                                                                                             |
| **R**                                            | 4.4.0                    | same / `r-base=4.4.0`                                                    | All steps                                                                                                                 | R language runtime                                                                                                                                                                                               | `R_MAX_VSIZE=48Gb` set in environment                                                                         |
| **Bioconductor**                                 | 3.19                     | (base image)                                                             | All steps                                                                                                                 | Bioconductor package ecosystem                                                                                                                                                                                   | —                                                                                                             |
| **nf-schema**                                    | 2.1.1                    | Nextflow plugin                                                          | Pipeline initialisation                                                                                                   | Parameter validation against JSON schema                                                                                                                                                                         | —                                                                                                             |

---

## SECTION 5 — PARAMETERS

| Parameter                   | Type                                 | Default       | Required | Description                                                                                  | Affects                                                   |
| --------------------------- | ------------------------------------ | ------------- | -------- | -------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `--input`                   | string (file-path, .tsv)             | null          | **Yes**  | Samplesheet TSV (Sample_Name, Sentrix_ID, Sentrix_Position)                                  | IMPORT_IDAT                                               |
| `--idat_dir`                | string (directory-path)              | null          | **Yes**  | Directory containing raw IDAT files                                                          | IMPORT_IDAT, SESAME_POOBAH_QC                             |
| `--meta_file`               | string (file-path, .csv/.tsv)        | null          | **Yes**  | Sample metadata CSV with clinical/phenotype info                                             | SEX_QC, ALIGN_META, COVARIATE_PCA, DMP_LIMMA, DMR_DMRCATE |
| `--outdir`                  | string (directory-path)              | `results`     | **Yes**  | Output directory                                                                             | All publishDir directives                                 |
| `--keep_groups`             | string (comma-separated)             | null          | **Yes**  | Groups to retain for EWAS; ≥2 required; 'CONTROL' is reference level                         | ALIGN_META, DMR_DMRCATE                                   |
| `--annotation`              | string                               | `20a1.hg38`   | No       | EPICv2 annotation manifest version                                                           | IMPORT_IDAT, FILTER_XY_NONCG_BEADS                        |
| `--genome_build`            | string (enum: hg38, hg19)            | `hg38`        | No       | Reference genome build for DMR coordinates                                                   | DMR_DMRCATE                                               |
| `--p_det`                   | number [0,1]                         | `0.05`        | No       | SeSAMe pOOBAH detection p-value threshold                                                    | SESAME_POOBAH_QC                                          |
| `--fail_sample_fr`          | number [0,1]                         | `0.10`        | No       | Max fraction of failed probes per sample before removal                                      | SESAME_POOBAH_QC                                          |
| `--fail_probe_fr`           | number [0,1]                         | `0.10`        | No       | Max fraction of failed samples per probe before removal                                      | SESAME_POOBAH_QC                                          |
| `--sesame_prep`             | string                               | `QCDPB`       | No       | SeSAMe preprocessing pipeline code (Q=quality mask, C=channel, D=dye bias, P=pOOBAH, B=NOOB) | SESAME_POOBAH_QC                                          |
| `--do_rmsnpandch`           | boolean                              | `true`        | No       | Remove SNP-overlapping and cross-hybridising probes                                          | SPLIT_COLLAPSE                                            |
| `--snp_dist_bp`             | integer ≥0                           | `2`           | No       | Distance in bp from SNP to flag probe                                                        | SPLIT_COLLAPSE                                            |
| `--snp_maf`                 | number [0,1]                         | `0.01`        | No       | Minor allele frequency threshold for SNP probe filtering                                     | SPLIT_COLLAPSE                                            |
| `--do_estimate_cellcomp`    | boolean                              | `true`        | No       | Estimate blood cell composition using FlowSorted.Blood.EPIC                                  | CELL_COMP; ALIGN_META (cell covariates)                   |
| `--filter_beads`            | boolean                              | `true`        | No       | Filter probes with insufficient bead counts                                                  | FILTER_XY_NONCG_BEADS                                     |
| `--bead_min`                | integer ≥1                           | `3`           | No       | Minimum beads required per probe per sample                                                  | FILTER_XY_NONCG_BEADS                                     |
| `--bead_fail_frac`          | number [0,1]                         | `0.05`        | No       | Max fraction of samples below bead_min for probe retention                                   | FILTER_XY_NONCG_BEADS                                     |
| `--filter_xy`               | boolean                              | `true`        | No       | Remove sex chromosome probes (chrX, chrY)                                                    | FILTER_XY_NONCG_BEADS                                     |
| `--filter_noncg`            | boolean                              | `true`        | No       | Remove non-CpG probes (retain only `cg*`)                                                    | FILTER_XY_NONCG_BEADS                                     |
| `--export_csv_matrices`     | boolean                              | `false`       | No       | Export final beta/M-value matrices as gzip CSV                                               | FINAL_EXPORTS, ALIGN_META                                 |
| `--save_final_matrices_rds` | boolean                              | `false`       | No       | Save final matrices as RDS to results                                                        | FINAL_EXPORTS                                             |
| `--matrix_rds_compression`  | string (enum: xz, gzip, bzip2, none) | `xz`          | No       | Compression for RDS matrix files                                                             | FINAL_EXPORTS                                             |
| `--save_intermediates`      | boolean                              | `false`       | No       | Save rg0.rds, beta_norm.rds, beta3.rds to results                                            | IMPORT_IDAT, SESAME_POOBAH_QC, FILTER_XY_NONCG_BEADS      |
| `--cofactors`               | string (comma-separated)             | `sex,age,bmi` | No       | Covariates for design matrix (sex, age, bmi, or any metadata column)                         | ALIGN_META, COVARIATE_PCA, DMP_LIMMA, DMR_DMRCATE         |
| `--epicv2_remap`            | boolean                              | `true`        | No       | Remap EPICv2 replicate probes for DMRcate                                                    | DMR_DMRCATE                                               |
| `--posreps_strategy`        | string (enum: mean, first, random)   | `mean`        | No       | Strategy for collapsing EPICv2 positional replicates                                         | DMR_DMRCATE                                               |
| `--n_probes_pca`            | integer ≥1000                        | `50000`       | No       | Number of most-variable probes for PCA                                                       | COVARIATE_PCA                                             |
| `--n_pcs`                   | integer ≥1                           | `10`          | No       | Number of PCs to test for covariate association                                              | COVARIATE_PCA                                             |
| `--run_covariate_pca`       | boolean                              | `true`        | No       | Run covariate PCA analysis                                                                   | COVARIATE_PCA                                             |
| `--run_champ_svd`           | boolean                              | `true`        | No       | Run ChAMP SVD analysis                                                                       | COVARIATE_PCA                                             |
| `--run_dmp`                 | boolean                              | `true`        | No       | Run DMP analysis with limma                                                                  | DMP_LIMMA                                                 |
| `--run_dmr`                 | boolean                              | `true`        | No       | Run DMR analysis with DMRcate                                                                | DMR_DMRCATE                                               |
| `--run_cometh`              | boolean                              | `true`        | No       | Run coMethDMR alongside DMRcate for consensus DMRs                                           | DMR_DMRCATE                                               |
| `--include_array`           | boolean                              | `false`       | No       | Include Sentrix_Position (Array) as batch covariate                                          | ALIGN_META                                                |
| `--publish_dir_mode`        | string                               | `copy`        | No       | Mode for saving results to output directory                                                  | All publishDir                                            |
| `--validate_params`         | boolean                              | `true`        | No       | Validate parameters against schema at runtime                                                | Pipeline init                                             |
| `--email`                   | string                               | null          | No       | Email address for completion summary                                                         | PIPELINE_COMPLETION                                       |
| `--email_on_fail`           | string                               | null          | No       | Email address for failure notification only                                                  | PIPELINE_COMPLETION                                       |
| `--hook_url`                | string                               | null          | No       | Webhook URL for messaging notifications                                                      | PIPELINE_COMPLETION                                       |
| `--monochrome_logs`         | boolean                              | `false`       | No       | Disable colour in log output                                                                 | All logging                                               |
| `--plaintext_email`         | boolean                              | `false`       | No       | Send plain-text instead of HTML email                                                        | PIPELINE_COMPLETION                                       |

---

## SECTION 6 — QUALITY CONTROL STRATEGY

### QC Steps and Metrics

**Step 1 — Raw fluorescence intensity QC (`RAW_INTENSITIES_QC`)**

- Method: `minfi::preprocessRaw()` + `minfi::getQC()`
- Metrics: Per-sample median methylated (mMed) and unmethylated (uMed) channel intensities
- Threshold: mMed < 10.5 OR uMed < 10.5 → sample flagged (informational; no removal at this stage)
- Output: Scatter plot, QC TSV, flagged sample list

**Step 2 — SeSAMe pOOBAH detection p-value QC (`SESAME_POOBAH_QC`)**

- Method: `sesame::openSesame()` with pOOBAH masking; probes failing detection p-value are set to NA
- Metrics: per-sample masking rate CSV, per-probe masking rate CSV, distribution histogram
- Sample threshold: fail rate > `--fail_sample_fr` (default 10%) → removed → `dropped_samples.csv`
- Probe threshold: fail rate > `--fail_probe_fr` (default 10%) → removed from beta matrix

**Step 3 — Control probe QC (`QC_REPORT_CONTROLS`)**

- Method: `minfi::qcReport()`, `minfi::controlStripPlot()`
- Control categories: STAINING, EXTENSION, HYBRIDIZATION, TARGET REMOVAL, SPECIFICITY I, SPECIFICITY II, NON-POLYMORPHIC, NEGATIVE, BISULFITE CONVERSION I, BISULFITE CONVERSION II
- No automatic pass/fail thresholds; outputs for manual visual inspection

**Step 4 — Sample identity QC (`SNP_HEATMAP`)**

- Method: Pairwise Pearson correlation of 59 SNP probe beta values; hierarchical clustering heatmap
- No automated flag threshold; visual inspection for unexpected similarity between distinct samples

**Step 5 — Beta-value distribution QC (`DENSITY_PLOTS`)**

- Method: Overlay density plots of raw and SeSAMe-normalized beta values
- Purpose: Identify outlier samples; assess normalization effect

**Step 6 — Sex concordance QC (`SEX_QC`)**

- Method: `minfi::preprocessNoob()` + `minfi::getSex()` using X/Y chromosome median intensities
- Metrics: Predicted sex (M/F), xMed, yMed
- Concordance: Predicted sex compared to `sex` column in metadata; discordant samples listed in `sex_discordant_samples.csv`; flagging only, no automatic removal

**Step 7 — Probe-level filters (`FILTER_XY_NONCG_BEADS`)**

- (a) Bead count: probes with < `--bead_min` (default 3) beads in > `--bead_fail_frac` (default 5%) of samples → removed
- (b) Sex chromosome removal: chrX, chrY probes → removed
- (c) Non-CpG removal: probes not starting with `cg` → removed
- Output: filtering summary TSV with probe counts at each stage

**Step 8 — DMP genomic inflation factor (`DMP_LIMMA`)**

- Metric: Lambda GC = `median(qchisq(1-p, df=1)) / qchisq(0.5, df=1)` per contrast
- Reported in `DMP_summary_counts_and_lambda.csv`

**MultiQC:** NOT USED (explicitly listed in `.nf-core.yml` `skip_features`).

---

## SECTION 7 — REFERENCE DATA AND GENOME HANDLING

**Array manifest and annotation:** Built into Bioconductor packages — no user-supplied reference files required.

- `IlluminaHumanMethylationEPICv2manifest` (v1.0.0): probe definitions, channel assignments, probe types (Type I/II)
- `IlluminaHumanMethylationEPICv2anno.20a1.hg38` (v1.0.0): genomic annotations for EPICv2 probes (GRCh38/hg38 default, hg19 legacy)
- `--annotation` (default `20a1.hg38`) and `--genome_build` (default `hg38`) control versions and coordinate system

**SeSAMe quality masks:** Pre-computed within `sesameData` package; `sesameDataCache()` ensures local availability.

**Cell composition reference:** `FlowSorted.Blood.EPIC::IDOLOptimizedCpGs` and `IDOLOptimizedCpGs.compTable` — built-in IDOL reference profiles for 6 blood cell types; no user file required.

**SNP filtering database:** Embedded in `DMRcate::rmSNPandCH()`.

**Required user-supplied files:**

1. IDAT files (red and green channel per sample)
2. Samplesheet TSV
3. Metadata CSV

---

## SECTION 8 — COMPUTATIONAL REQUIREMENTS

### Resource labels (`conf/base.config`)

| Label                 | CPUs         | Memory           | Time           |
| --------------------- | ------------ | ---------------- | -------------- |
| `process_single`      | 1            | 6 GB × attempt   | 4 h × attempt  |
| `process_low`         | 2 × attempt  | 12 GB × attempt  | 4 h × attempt  |
| `process_medium`      | 6 × attempt  | 36 GB × attempt  | 8 h × attempt  |
| `process_high`        | 12 × attempt | 72 GB × attempt  | 16 h × attempt |
| `process_long`        | (inherited)  | (inherited)      | 20 h × attempt |
| `process_high_memory` | (inherited)  | 200 GB × attempt | (inherited)    |

### Per-process resource overrides

| Process               | CPUs | Memory          | Time          |
| --------------------- | ---- | --------------- | ------------- |
| IMPORT_IDAT           | 2    | 8 GB × attempt  | 4 h × attempt |
| RAW_INTENSITIES_QC    | 2    | 8 GB × attempt  | 2 h × attempt |
| LOG_VERSIONS          | 1    | 4 GB × attempt  | 1 h × attempt |
| SESAME_POOBAH_QC      | 4    | 16 GB × attempt | 8 h × attempt |
| QC_REPORT_CONTROLS    | 2    | 8 GB × attempt  | 4 h × attempt |
| SNP_HEATMAP           | 2    | 8 GB × attempt  | 2 h × attempt |
| DENSITY_PLOTS         | 2    | 8 GB × attempt  | 2 h × attempt |
| SEX_QC                | 2    | 12 GB × attempt | 4 h × attempt |
| FILTER_XY_NONCG_BEADS | 2    | 16 GB × attempt | 6 h × attempt |
| CELL_COMP             | 2    | 16 GB × attempt | 6 h × attempt |
| SPLIT_COLLAPSE        | 2    | 12 GB × attempt | 4 h × attempt |
| FINAL_EXPORTS         | 1    | 8 GB × attempt  | 2 h × attempt |
| ALIGN_META            | 1    | 8 GB × attempt  | 2 h × attempt |
| COVARIATE_PCA         | 2    | 8 GB × attempt  | 4 h × attempt |
| DMP_LIMMA             | 2    | 16 GB × attempt | 8 h × attempt |
| DMR_DMRCATE           | 2    | 16 GB × attempt | 8 h × attempt |

**Global hard caps:** 16 CPUs, 128 GB memory, 240 hours

**R memory:** `R_MAX_VSIZE=48Gb` set via `env {}` block in `nextflow.config`

**Error strategy:** `errorStrategy = { task.exitStatus in ((130..145) + 104 + 175) ? 'retry' : 'finish' }`; `maxRetries = 1`; resources scale linearly (2× on retry).

**Supported executors:** Any Nextflow-compatible executor (local, SLURM, SGE, PBS, LSF, AWS Batch, Google Cloud Batch, etc.) via institutional profiles.

**Available config profiles:** `docker`, `singularity`, `podman`, `shifter`, `charliecloud`, `apptainer`, `conda`, `mamba`, `wave`, `arm64`, `emulate_amd64`, `gpu`, `debug`, `test`, `test_full`

---

## SECTION 9 — CONTAINERIZATION AND REPRODUCIBILITY

**Supported container technologies:** Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer; also Conda/Mamba.

**Container registry:** `quay.io` (default for all container engines in `nextflow.config`)

**Pipeline container:** `quay.io/nf-core/methylarray:1.0.0dev` — used by all 16 local R modules.

### Per-module conda packages (pinned versions)

| Module                | Key conda packages                                                                                                                                                                                               |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| IMPORT_IDAT           | `bioconductor-minfi=1.48.0`, `bioconductor-illuminahumanmethylationepicv2manifest=1.0.0`, `bioconductor-illuminahumanmethylationepicv2anno.20a1.hg38=1.0.0`, `r-readr=2.1.5`, `r-optparse=1.7.5`, `r-base=4.4.0` |
| RAW_INTENSITIES_QC    | `bioconductor-minfi=1.48.0`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                                                  |
| SESAME_POOBAH_QC      | `bioconductor-sesame=1.18.0`, `bioconductor-sesamedata=1.18.0`, `bioconductor-minfi=1.48.0`, `r-readr=2.1.5`                                                                                                     |
| QC_REPORT_CONTROLS    | `bioconductor-minfi=1.48.0`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                                                  |
| SNP_HEATMAP           | `bioconductor-minfi=1.48.0`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                                                  |
| DENSITY_PLOTS         | `bioconductor-minfi=1.48.0`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                                                  |
| SEX_QC                | `bioconductor-minfi=1.48.0`, `r-readr=2.1.5`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                                 |
| FILTER_XY_NONCG_BEADS | `bioconductor-minfi=1.48.0`, `bioconductor-illuminahumanmethylationepicv2anno.20a1.hg38=1.0.0`, `r-dplyr=1.1.4`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                              |
| CELL_COMP             | `bioconductor-minfi=1.48.0`, `bioconductor-flowsorted.blood.epic=2.4.0`, `r-readr=2.1.5`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                     |
| SPLIT_COLLAPSE        | `bioconductor-dmrcate=2.16.0`, `bioconductor-sesame=1.18.0`, `bioconductor-illuminahumanmethylationepicv2anno.20a1.hg38=1.0.0`, `r-optparse=1.7.5`, `r-base=4.4.0`                                               |
| FINAL_EXPORTS         | `r-readr=2.1.5`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                                                              |
| ALIGN_META            | `r-readr=2.1.5`, `r-dplyr=1.1.4`, `r-tibble=3.2.1`, `r-optparse=1.7.5`, `r-base=4.4.0`                                                                                                                           |
| COVARIATE_PCA         | `bioconductor-champ=2.32.0`, `bioconductor-sva=3.50.0`, `r-readr=2.1.5`, `r-dplyr=1.1.4`, `r-matrixstats=1.3.0`, `r-ggplot2=3.5.1`, `r-optparse=1.7.5`, `r-base=4.4.0`                                           |
| DMP_LIMMA             | `bioconductor-limma=3.58.0`, `bioconductor-missmethyl=1.36.0`, `bioconductor-sva=3.50.0`, `r-ggplot2=3.5.1`, `r-dplyr=1.1.4`, `r-readr=2.1.5`, `r-tibble=3.2.1`, `r-optparse=1.7.5`, `r-base=4.4.0`              |
| DMR_DMRCATE           | `bioconductor-dmrcate=2.16.0`, `bioconductor-minfi=1.48.0`, `bioconductor-limma=3.58.0`, `bioconductor-sva=3.50.0`, `r-readr=2.1.5`, `r-ggplot2=3.5.1`, `r-dplyr=1.1.4`, `r-optparse=1.7.5`, `r-base=4.4.0`      |
| LOG_VERSIONS          | `r-base=4.4.0`, `r-optparse=1.7.5`                                                                                                                                                                               |

**Conda channels:** `conda-forge`, `bioconda` (no `defaults` channel)

**Software versions:** Pinned exactly in all `environment.yml` files; no build hashes (nf-core compliant).

**Reproducibility measures:**

- Nextflow `-resume` for task-level caching
- Containers ensure identical software environments
- `set -euo pipefail` shell options enforced
- `R_PROFILE_USER` and `R_ENVIRON_USER` set to empty paths (no local R config interference)
- PCA uses `set.seed(1)` for deterministic output

**Test profile:** `conf/test.config` — 4 CPUs, 15 GB, 1 h; slow optional steps disabled; minimal dataset from `nf-core/test-datasets/methylarray/`.

**Software version logging:** `LOG_VERSIONS` module calls `log_versions.R` which queries `packageVersion()` for 18 packages; all modules emit `versions.yml` collated into `nf_core_methylarray_software_versions.yml`.

---

## SECTION 10 — CITATIONS AND REFERENCES

### nf-core framework

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. _Nat Biotechnol._ 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PMID: 32055031.

### Nextflow

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. _Nat Biotechnol._ 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PMID: 28398311.

### minfi

> Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. _Bioinformatics._ 2014 May 15;30(10):1363-9. doi: 10.1093/bioinformatics/btu049. PMID: 24478339.

### SeSAMe

> Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. _Nucleic Acids Res._ 2018 Nov 2;46(20):e123. doi: 10.1093/nar/gky691. PMID: 30085087.

### limma

> Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. _Nucleic Acids Res._ 2015 Apr 20;43(7):e47. doi: 10.1093/nar/gkv007. PMID: 25605792.

### DMRcate

> Peters TJ, Buckley MJ, Statham AL, Pidsley R, Samaras K, R VL, Clark SJ, Molloy PL. De novo identification of differentially methylated regions in the human genome. _Epigenetics Chromatin._ 2015 Jan 27;8:6. doi: 10.1186/1756-8935-8-6. PMID: 25972926.

### ChAMP

> Tian Y, Morris TJ, Webster AP, Yang Z, Beck S, Feber A, Teschendorff AE. ChAMP: updated methylation analysis pipeline for Illumina BeadChips. _Bioinformatics._ 2017 Dec 15;33(24):3982-3984. doi: 10.1093/bioinformatics/btx513. PMID: 28961746.

### FlowSorted.Blood.EPIC

> Salas LA, Zhang Z, Koestler DC, Butler RA, Hansen KD, Molinaro AM, Wiencke JK, Kelsey KT, Christensen BC. Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling. _Nat Commun._ 2022 Feb 3;13(1):761. doi: 10.1038/s41467-021-27864-7. PMID: 35115528.

### Houseman deconvolution

> Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, Nelson HH, Wiencke JK, Kelsey KT. DNA methylation arrays as surrogate measures of cell mixture distribution. _BMC Bioinformatics._ 2012 Mar 28;13:86. doi: 10.1186/1471-2105-13-86. PMID: 22568884.

### ComBat / sva

> Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. _Biostatistics._ 2007 Jan;8(1):118-27. doi: 10.1093/biostatistics/kxj037. PMID: 16632515.

### R

> R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.

### Software packaging and containers

- **Bioconda:** Grüning B et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. _Nat Methods._ 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PMID: 29967506.
- **BioContainers:** da Veiga Leprevost F et al. BioContainers: an open-source and community-driven framework for software standardization. _Bioinformatics._ 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PMID: 28379341.
- **Docker:** Merkel D. Docker: lightweight linux containers for consistent development and deployment. _Linux Journal._ 2014(239):2. doi: 10.5555/2600239.2600241.
- **Singularity:** Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. _PLoS One._ 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. PMID: 28494014.

**Pipeline's own citation:** [NOT FOUND — Zenodo DOI placeholder (`10.5281/zenodo.XXXXXXX`) only; no real DOI published yet]

---

## SECTION 11 — STATISTICAL AND ALGORITHMIC METHODS

### 1. Raw intensity QC (`raw_intensities_qc.R`)

- **Method:** `minfi::preprocessRaw()` + `minfi::getQC()` → per-sample median log2-intensities for methylated and unmethylated channels
- **Threshold:** Hardcoded at 10.5 (`badSampleCutoff=10.5` from minfi default); flagging only
- **Assumption:** High-quality samples have median intensities substantially above background

### 2. SeSAMe normalization and pOOBAH QC (`sesame_poobah_qc.R`)

- **Method:** `sesame::openSesame()` with preparation code `QCDPB`:
  - **Q** — Genomic quality masking: pre-computed masks from `sesameData` (bead count < 3, probes in common genomic deletions) → NA
  - **C** — Colour channel inference: corrects Infinium I probe colour channel switching artefacts
  - **D** — Non-linear dye bias correction
  - **P** — pOOBAH p-value masking: out-of-band (OOB) signal as null distribution; probes with p > `--p_det` (default 0.05) → NA
  - **B** — NOOB background correction: negative control probe-based subtraction (Normal-Exponential Out-Of-Band, Triche 2013 method)
- **No probe-level FDR correction** — raw pOOBAH p-value used as per-observation quality mask
- **Assumptions:** OOB signal is a valid null distribution; background fluorescence is modelled as exponential; dye imbalance is systematic across probes

### 3. Sex prediction (`sex_qc.R`)

- **Method:** `minfi::preprocessNoob()` → `minfi::mapToGenome()` → `minfi::getSex()`: 2D threshold on median X intensity (xMed) vs median Y intensity (yMed)
- **Concordance:** Fuzzy-matched against metadata `sex` column (M/F/Male/Female/0/1 variants accepted)
- **Assumptions:** X/Y chromosome probe intensities reliably predict biological sex; XY males show high yMed, XX females show low yMed/high xMed

### 4. Bead count filtering (`filter_xy_noncg_beads.R`)

- **Method:** `minfi::getNBeads()` → address-to-probe-name mapping via `minfi::getProbeInfo()`; Type I probes: min(AddressA, AddressB); Type II probes: AddressA
- **Threshold:** `mean(bead_count < bead_min) > bead_fail_frac` → probe removed (default: bead_min=3, bead_fail_frac=0.05)
- **Assumption:** < 3 beads per probe → unreliable signal due to insufficient physical measurement redundancy

### 5. Cell composition deconvolution (`cell_comp.R`)

- **Method:** Houseman constraint projection (`FlowSorted.Blood.EPIC::projectCellType_CP`) using IDOL optimized CpG set
- **Reference:** 450 IDOL CpGs (`IDOLOptimizedCpGs`); reference composition matrix (`IDOLOptimizedCpGs.compTable`) for 6 cell types (CD8T, CD4T, NK, Bcell, Mono, Neu, Gran)
- **Input processing:** EPICv2 positional replicates collapsed by prefix (mean) before IDOL intersection; `preprocessNoob` normalization applied first
- **Constraints:** Non-negativity only (`nonnegative=TRUE`, `lessThanOne=FALSE`)
- **Minimum requirement:** ≥200 IDOL CpG overlap in input data
- **Assumptions:** Blood cell methylation is adequately captured by IDOL reference profiles derived from EPIC arrays; profiles are transferable to EPICv2

### 6. SNP and cross-hybridisation probe filtering (`split_collapse.R`)

- **Method:** `DMRcate::rmSNPandCH()` using internal SNP database
- **SNP threshold:** Probes within `--snp_dist_bp` (default 2 bp) of a SNP with MAF ≥ `--snp_maf` (default 0.01) → removed
- **Cross-hybridisation:** Removed for DMP branch only; DMR branch retains these probes
- **EPICv2 replicate collapsing:** `sesame::betasCollapseToPfx()` with `posreps_strategy` (mean/first/random) for DMP branch

### 7. Covariate PCA (`covariate_pca.R`)

- **Method:** `prcomp(t(beta_sub), center=TRUE, scale.=FALSE)` on top `n_probes_pca` (default 50,000) most variable CpGs selected by `matrixStats::rowVars()`
- **Missing value imputation:** Row means used to impute NA values before PCA
- **Association testing:** `lm(PC ~ covariate)` — F-test p-values for continuous covariates; one-way ANOVA for categorical; results as -log10(p) matrix
- **Covariate balance tests:** Chi-square (`chisq.test`) for Plate × Group; one-way ANOVA for Age, BMI, cell fractions by group
- **Seed:** `set.seed(1)`
- **ChAMP SVD:** `ChAMP::champ.SVD()` applied separately on beta values with ≤1% NA rate

### 8. DMP analysis — limma (`dmp_limma.R`)

- **Response variable:** M-values: `M = log2(β_clipped / (1 - β_clipped))` with clipping ε = 1e-6
- **Design matrix:** Cell-means coding (`~ 0 + group + cofactors`); all pairwise contrasts generated automatically via `combn(groups, 2)` + `makeContrasts`
- **`eBayes` settings:** `trend=TRUE` (accounts for mean-variance relationship in M-values), `robust=TRUE` (reduces influence of outlier probes)
- **NA handling:** Probes with > 1% non-finite M-values removed; remaining NAs → row-mean imputation
- **Multiple testing correction:** Benjamini-Hochberg (BH) FDR via `limma::topTable(adjust.method="BH")`
- **Delta-beta:** Arithmetic mean of beta values per group; difference added to DMP table
- **Significance summary:** n DMPs at FDR < 0.05 and lambda_GC per contrast
- **Three analysis modes run in parallel:**
  - `norm_only`: baseline model (group + cofactors + cell composition)
  - `plate_in_model`: adds `Plate` (Sentrix_ID) as a fixed covariate
  - `combat`: ComBat batch correction (`sva::ComBat`, `par.prior=TRUE`) on M-values using Plate as batch; then SVA (`sva::sva`) with "be" and "leek" k-selection; backoff loop ensures full-rank design
- **Cell composition covariate:** All estimated cell types included as covariates minus one (Neu dropped to avoid collinearity)
- **Assumptions:** Linear additive model for M-values; homoscedastic variance after eBayes shrinkage; ComBat assumes Gaussian additive batch effects

### 9. DMR analysis — DMRcate (`dmr_dmrcate.R`)

- **Step 1 — CpG annotation:** `DMRcate::cpg.annotate()`: `analysis.type="differential"`, `what="M"`, `arraytype="EPICv2"`, `epicv2Filter=posreps_strategy`, `epicv2Remap=TRUE`, `fdr=0.05`; uses limma-derived statistics from same design
- **Step 2 — DMR calling:** `DMRcate::dmrcate(ann, lambda=1000, C=2)`: Gaussian kernel smoother (bandwidth 1000 bp); minimum 2 CpGs per DMR
- **Step 3 — Range extraction:** `DMRcate::extractRanges()`: genomic coordinates, CpG count, min FDR, Stouffer combined p-value, mean methylation difference
- **Input:** Uncollapsed EPICv2 beta matrix (`beta_for_DMR_uncollapsed.rds`); M-values computed internally with same logit transform
- **Three analysis modes:** Same as DMP (norm_only, plate_in_model, combat with SVA)
- **Minimum requirement:** ≥6 samples before DMR analysis

### 10. coMethDMR analysis (optional) (`dmr_dmrcate.R`)

- **Method:** `coMethDMR::CoMethAllRegions()` on each pairwise contrast; binary phenotype coding (1=case, 0=control); serial execution (`BiocParallel::SerialParam()`)
- **FDR threshold:** < 0.05 to retain coMethDMR regions
- **Consensus:** `GenomicRanges::findOverlaps(minoverlap=1)` between DMRcate and coMethDMR results; DMRcate statistics reported for overlapping DMRs

### Summary of multiple testing correction methods

| Step                     | Method                                               | Software         |
| ------------------------ | ---------------------------------------------------- | ---------------- |
| pOOBAH probe masking     | No correction; threshold p < 0.05 per observation    | SeSAMe           |
| DMP analysis             | Benjamini-Hochberg (BH) FDR across all tested probes | limma `topTable` |
| DMR analysis             | FDR via DMRcate's Stouffer method per region         | DMRcate          |
| coMethDMR                | FDR per region (coMethDMR internal)                  | coMethDMR        |
| PC-covariate association | Raw p-values reported as -log10; no FDR applied      | prcomp + lm      |

---

## SECTION 12 — OUTPUT FILES

### Primary results

| File / Directory                                                          | Type      | Content                                                                          | Producing step |
| ------------------------------------------------------------------------- | --------- | -------------------------------------------------------------------------------- | -------------- |
| `final_exports/bVals_final.rds`                                           | R matrix  | Final normalized beta values (probes × samples)                                  | FINAL_EXPORTS  |
| `final_exports/mVals_final.rds`                                           | R matrix  | Final M-values (logit of beta)                                                   | FINAL_EXPORTS  |
| `final_exports/probes_final.csv`                                          | CSV       | Probe manifest (probe IDs retained in final matrices)                            | FINAL_EXPORTS  |
| `final_exports/samples_final.csv`                                         | CSV       | Sample manifest                                                                  | FINAL_EXPORTS  |
| `final_exports/qc_summary.txt`                                            | text      | Final probe count × sample count summary                                         | FINAL_EXPORTS  |
| `align_meta/meta_aligned.csv`                                             | CSV       | Aligned metadata with all covariates, predicted sex, cell composition            | ALIGN_META     |
| `align_meta/bVals_aligned.rds`                                            | R matrix  | Beta matrix subset to aligned metadata samples                                   | ALIGN_META     |
| `align_meta/mVals_aligned.rds`                                            | R matrix  | M-value matrix subset to aligned metadata samples                                | ALIGN_META     |
| `cell_comp/cell_counts_blood.csv`                                         | CSV       | Estimated cell type proportions per sample                                       | CELL_COMP      |
| `covariate_pca/covariates_PC_association_minuslog10p.csv`                 | CSV       | -log10(p) matrix: covariates × PCs                                               | COVARIATE_PCA  |
| `covariate_pca/covariate_checks.txt`                                      | text      | Covariate balance statistics (chi-square, ANOVA)                                 | COVARIATE_PCA  |
| `covariate_pca/ChAMP_SVD/`                                                | directory | ChAMP SVD analysis outputs                                                       | COVARIATE_PCA  |
| `dmp_limma/DMP_limma__norm_only/DMP_<A>_vs_<B>.csv`                       | CSV       | DMP results per contrast: logFC, t, P.Value, adj.P.Val, B, deltaBeta             | DMP_LIMMA      |
| `dmp_limma/DMP_limma__plate_in_model/DMP_<A>_vs_<B>.csv`                  | CSV       | DMP results with Plate as fixed covariate                                        | DMP_LIMMA      |
| `dmp_limma/DMP_limma__combat/DMP_<A>_vs_<B>.csv`                          | CSV       | DMP results after ComBat + SVA batch correction                                  | DMP_LIMMA      |
| `dmp_limma/*/DMP_summary_counts_and_lambda.csv`                           | CSV       | n DMPs at FDR < 0.05 and lambda_GC per contrast                                  | DMP_LIMMA      |
| `dmr_dmrcate/DMR_DMRcate_EPICv2__*/DMRcate_<A>_vs_<B>_ranges.csv`         | CSV       | DMR results: seqnames, start, end, no.cpgs, min_smoothed_fdr, meandiff, Stouffer | DMR_DMRCATE    |
| `dmr_dmrcate/DMR_DMRcate_EPICv2__*/DMRcate_summary_nDMR_per_contrast.csv` | CSV       | Count of significant DMRs per contrast                                           | DMR_DMRCATE    |
| `dmr_dmrcate/DMR_coMethDMR__*/coMethDMR_<A>_vs_<B>_ranges.csv`            | CSV       | coMethDMR results (if `--run_cometh true`)                                       | DMR_DMRCATE    |
| `dmr_dmrcate/DMR_consensus__*/consensus_<A>_vs_<B>_ranges.csv`            | CSV       | Consensus DMRs from DMRcate ∩ coMethDMR overlap                                  | DMR_DMRCATE    |

### QC outputs

| File                                                       | Type     | Content                                                 | Step                  |
| ---------------------------------------------------------- | -------- | ------------------------------------------------------- | --------------------- |
| `raw_intensities_qc/01_plotQC_raw_intensities.png`         | PNG      | Scatter plot mMed vs uMed per sample                    | RAW_INTENSITIES_QC    |
| `raw_intensities_qc/qc_raw.tsv`                            | TSV      | Per-sample mMed and uMed values                         | RAW_INTENSITIES_QC    |
| `raw_intensities_qc/low_intensity_samples.tsv`             | TSV      | Flagged low-intensity samples                           | RAW_INTENSITIES_QC    |
| `sesame_poobah_qc/02_hist_failrate_per_sample_masking.png` | PNG      | Histogram of per-sample pOOBAH masking rates            | SESAME_POOBAH_QC      |
| `sesame_poobah_qc/masking_failrate_per_sample.csv`         | CSV      | Per-sample pOOBAH fail rate                             | SESAME_POOBAH_QC      |
| `sesame_poobah_qc/masking_failrate_per_probe.csv`          | CSV      | Per-probe pOOBAH fail rate                              | SESAME_POOBAH_QC      |
| `sesame_poobah_qc/dropped_samples.csv`                     | CSV      | Samples removed by pOOBAH QC                            | SESAME_POOBAH_QC      |
| `sesame_poobah_qc/targets_filtered.tsv`                    | TSV      | Targets table after sample removal                      | SESAME_POOBAH_QC      |
| `qc_report_controls/00_qcReport.pdf`                       | PDF      | Comprehensive minfi QC report                           | QC_REPORT_CONTROLS    |
| `qc_report_controls/00_controlStrip_*.png`                 | PNG (×8) | Strip plots per Illumina control type                   | QC_REPORT_CONTROLS    |
| `qc_report_controls/03_control_bisulfite_I.png`            | PNG      | Type I bisulfite conversion efficiency                  | QC_REPORT_CONTROLS    |
| `qc_report_controls/04_control_bisulfite_II.png`           | PNG      | Type II bisulfite conversion efficiency                 | QC_REPORT_CONTROLS    |
| `snp_heatmap/05_heatmap_snp_correlation.png`               | PNG      | SNP probe correlation heatmap                           | SNP_HEATMAP           |
| `density_plots/06_density_beta_raw.png`                    | PNG      | Raw beta density overlay                                | DENSITY_PLOTS         |
| `density_plots/07_density_beta_sesame.png`                 | PNG      | SeSAMe-normalized beta density overlay                  | DENSITY_PLOTS         |
| `sex_qc/08_plotSex_predicted.png`                          | PNG      | Predicted sex scatter plot (X vs Y intensity)           | SEX_QC                |
| `sex_qc/qc_predicted_sex.csv`                              | CSV      | Predicted sex, X/Y intensities, reported sex per sample | SEX_QC                |
| `sex_qc/sex_discordant_samples.csv`                        | CSV      | Samples with sex prediction/metadata discordance        | SEX_QC                |
| `filter_xy_noncg_beads/filtering_summary.tsv`              | TSV      | Probe counts at each filtering stage                    | FILTER_XY_NONCG_BEADS |

### Intermediate outputs (optional, default: not saved)

| File                               | Condition                    | Content                             |
| ---------------------------------- | ---------------------------- | ----------------------------------- |
| `import_idat/rg0.rds`              | `--save_intermediates true`  | Raw RGChannelSet                    |
| `sesame_poobah_qc/beta_norm.rds`   | `--save_intermediates true`  | Post-SeSAMe, pre-filter beta matrix |
| `filter_xy_noncg_beads/beta3.rds`  | `--save_intermediates true`  | Post-filter beta matrix             |
| `final_exports/bVals_final.csv.gz` | `--export_csv_matrices true` | Beta matrix as gzip CSV             |
| `final_exports/mVals_final.csv.gz` | `--export_csv_matrices true` | M-value matrix as gzip CSV          |

### Pipeline execution info

| File                                                      | Content                                     |
| --------------------------------------------------------- | ------------------------------------------- |
| `pipeline_info/execution_timeline_*.html`                 | Per-task timeline visualization             |
| `pipeline_info/execution_report_*.html`                   | Execution summary report                    |
| `pipeline_info/execution_trace_*.txt`                     | Per-task resource usage trace               |
| `pipeline_info/pipeline_dag_*.html`                       | Workflow DAG visualization                  |
| `pipeline_info/nf_core_methylarray_software_versions.yml` | Collated software versions from all modules |
| `pipeline_info/package_versions.csv`                      | R package versions                          |
| `pipeline_info/sessionInfo.txt`                           | Full R session info                         |

---

## SECTION 13 — VALIDATION AND TEST DATA

**Test dataset:** Hosted on `nf-core/test-datasets` GitHub repository, branch `methylarray`.

**Test dataset URLs (`conf/test.config`):**

- Samplesheet: `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/methylarray/methylarray/samplesheet/samplesheet_test.tsv`
- Metadata: `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/methylarray/methylarray/metadata/meta_test.csv`
- IDAT files: same base URL, `idat/` subdirectory

**Full test dataset URLs (`conf/test_full.config`):**

- Samplesheet: `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/methylarray/methylarray/samplesheet/samplesheet_full.tsv`
- Metadata: `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/methylarray/methylarray/metadata/meta_full.csv`

**Test groups:** CONTROL, CCA (cholangiocarcinoma), PSC (primary sclerosing cholangitis), GBC (gallbladder cancer)

**Accession numbers:** [NOT FOUND — test data not linked to public accession; appears to be from an unpublished cohort stored on nf-core/test-datasets]

**Minimal test configuration:**

- Resource limits: 4 CPUs, 15 GB RAM, 1 hour
- Disabled: DMP analysis, DMR analysis, ChAMP SVD, covariate PCA, cell composition estimation
- Purpose: Verify end-to-end pipeline execution through FINAL_EXPORTS/ALIGN_META

**Benchmark validation:** [NOT FOUND — no published benchmark comparison documented]

---

## SECTION 14 — METHODS PARAGRAPH STUB

> All analyses were performed using the **nf-core/methylarray** pipeline (v1.0.0dev) implemented in Nextflow (≥25.04.0) following nf-core best practices [Ewels et al., *Nat Biotechnol*, 2020]. This pipeline processes Illumina EPICv2 DNA methylation array IDAT data through comprehensive quality control, normalization, and differential methylation analysis to identify differentially methylated positions (DMPs) and regions (DMRs) between disease groups. Briefly, raw IDAT files were imported and quality-assessed using minfi (v1.48.0) [Aryee et al., *Bioinformatics*, 2014]; samples with low fluorescence intensity (median log2-intensity < 10.5 in either channel) were flagged. Signal preprocessing and normalization were performed using SeSAMe (v1.18.0) [Zhou et al., *Nucleic Acids Res*, 2018] with the QCDPB pipeline comprising genomic quality masking, channel inference, non-linear dye bias correction, pOOBAH detection p-value masking (threshold p < 0.05), and NOOB background subtraction; samples with > 10% masked probes and probes failing in > 10% of samples were removed. Illumina internal control probe performance, bisulfite conversion efficiency, SNP-based sample identity, and predicted-vs-reported sex concordance were assessed. Probes on sex chromosomes and non-CpG probes were excluded, and probes with insufficient bead support (< 3 beads in > 5% of samples) were removed. Blood cell-type composition (CD8T, CD4T, NK, Bcell, Mono, Neu, Gran) was estimated by reference-based Houseman deconvolution using the IDOL optimized CpG set from FlowSorted.Blood.EPIC (v2.4.0) [Salas et al., *Nat Commun*, 2022; Houseman et al., *BMC Bioinformatics*, 2012]. SNP-overlapping and cross-hybridising probes were removed using DMRcate::rmSNPandCH (MAF ≥ 0.01, within 2 bp), and EPICv2 positional replicate probes were collapsed by mean. DMP analysis was performed on M-values using limma (v3.58.0) [Ritchie et al., *Nucleic Acids Res*, 2015] with empirical Bayes moderation (trend=TRUE, robust=TRUE); all pairwise group contrasts were tested with cell-type proportions, sex, age, and BMI as covariates; three parallel analyses were performed: baseline, plate-corrected, and batch-corrected using ComBat [Johnson et al., *Biostatistics*, 2007] with surrogate variable analysis. DMR analysis was performed using DMRcate (v2.16.0) [Peters et al., *Epigenetics Chromatin*, 2015] with EPICv2 remapping and Gaussian kernel smoothing (lambda = 1000 bp, minimum 2 CpGs per DMR). P-values were adjusted for multiple testing using the Benjamini-Hochberg method; significant DMPs were defined at FDR < 0.05. All software was run in containers using Docker/Singularity (image: `quay.io/nf-core/methylarray:1.0.0dev`). Full pipeline parameters are available in Supplementary Table X.

---

## RAW FILE INVENTORY

| File path                                                      | Status                                                    |
| -------------------------------------------------------------- | --------------------------------------------------------- |
| `main.nf`                                                      | Examined                                                  |
| `nextflow.config`                                              | Examined                                                  |
| `nextflow_schema.json`                                         | Examined                                                  |
| `README.md`                                                    | Examined                                                  |
| `CHANGELOG.md`                                                 | Examined                                                  |
| `CITATIONS.md`                                                 | Examined                                                  |
| `LICENSE`                                                      | Examined                                                  |
| `conf/base.config`                                             | Examined                                                  |
| `conf/modules.config`                                          | Examined                                                  |
| `conf/test.config`                                             | Examined                                                  |
| `conf/test_full.config`                                        | Examined                                                  |
| `workflows/methylarray.nf`                                     | Examined                                                  |
| `subworkflows/local/utils_nfcore_methylarray_pipeline/main.nf` | Examined                                                  |
| `assets/schema_input.json`                                     | Examined                                                  |
| `docs/usage.md`                                                | Examined                                                  |
| `docs/output.md`                                               | Examined                                                  |
| `modules.json`                                                 | Examined                                                  |
| `.nf-core.yml`                                                 | Examined                                                  |
| `environment.yml` (root)                                       | NOT FOUND — per-module environment.yml files used instead |
| `modules/local/import_idat/main.nf`                            | Examined                                                  |
| `modules/local/raw_intensities_qc/main.nf`                     | Examined                                                  |
| `modules/local/sesame_poobah_qc/main.nf`                       | Examined                                                  |
| `modules/local/qc_report_controls/main.nf`                     | Examined                                                  |
| `modules/local/snp_heatmap/main.nf`                            | Examined                                                  |
| `modules/local/density_plots/main.nf`                          | Examined                                                  |
| `modules/local/sex_qc/main.nf`                                 | Examined                                                  |
| `modules/local/filter_xy_noncg_beads/main.nf`                  | Examined                                                  |
| `modules/local/cell_comp/main.nf`                              | Examined                                                  |
| `modules/local/split_collapse/main.nf`                         | Examined                                                  |
| `modules/local/final_exports/main.nf`                          | Examined                                                  |
| `modules/local/align_meta/main.nf`                             | Examined                                                  |
| `modules/local/covariate_pca/main.nf`                          | Examined                                                  |
| `modules/local/dmp_limma/main.nf`                              | Examined                                                  |
| `modules/local/dmr_dmrcate/main.nf`                            | Examined                                                  |
| `modules/local/log_versions/main.nf`                           | Examined                                                  |
| `modules/local/import_idat/environment.yml`                    | Examined                                                  |
| `modules/local/sesame_poobah_qc/environment.yml`               | Examined                                                  |
| `modules/local/cell_comp/environment.yml`                      | Examined                                                  |
| `modules/local/split_collapse/environment.yml`                 | Examined                                                  |
| `modules/local/filter_xy_noncg_beads/environment.yml`          | Examined                                                  |
| `modules/local/covariate_pca/environment.yml`                  | Examined                                                  |
| `modules/local/dmp_limma/environment.yml`                      | Examined                                                  |
| `modules/local/dmr_dmrcate/environment.yml`                    | Examined                                                  |
| `bin/import_idat.R`                                            | Examined                                                  |
| `bin/raw_intensities_qc.R`                                     | Examined                                                  |
| `bin/sesame_poobah_qc.R`                                       | Examined                                                  |
| `bin/qc_report_controls.R`                                     | Examined                                                  |
| `bin/snp_heatmap.R`                                            | Examined                                                  |
| `bin/density_plots.R`                                          | Examined                                                  |
| `bin/sex_qc.R`                                                 | Examined                                                  |
| `bin/filter_xy_noncg_beads.R`                                  | Examined                                                  |
| `bin/cell_comp.R`                                              | Examined                                                  |
| `bin/split_collapse.R`                                         | Examined                                                  |
| `bin/final_exports.R`                                          | Examined                                                  |
| `bin/align_meta.R`                                             | Examined                                                  |
| `bin/covariate_pca.R`                                          | Examined                                                  |
| `bin/dmp_limma.R`                                              | Examined                                                  |
| `bin/dmr_dmrcate.R`                                            | Examined                                                  |
| `bin/log_versions.R`                                           | Examined                                                  |
| `bin/utils.R`                                                  | Examined                                                  |

---

## NOTES ON GAPS AND CAVEATS

- **Pipeline DOI:** Zenodo placeholder only (`10.5281/zenodo.XXXXXXX`); no real DOI assigned yet — update before publication.
- **Test data accession:** Not linked to any public repository (GEO, ArrayExpress, etc.); appears to be from an unpublished liver disease cohort.
- **coMethDMR version:** Not pinned in any `environment.yml`; version undetermined — pin before publication.
- **missMethyl:** Installed in DMP_LIMMA environment (v1.36.0) but not actively called in `dmp_limma.R`; reserved for future gene ontology/pathway enrichment analysis.
- **Container:** `quay.io/nf-core/methylarray:1.0.0dev` is a placeholder — build and publish real container before release.
- **ORCID identifiers:** Blank for both authors — add before submission.

# nf-core/methylarray: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [2026-03-15]

Initial release of nf-core/methylarray, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Full pipeline for Illumina EPIC v2 methylation array analysis (IDAT input).
- IMPORT_IDAT: Import raw IDAT files and build RGChannelSet + targets using minfi.
- RAW_INTENSITIES_QC: Quality control of raw signal intensities with density plots and low-intensity sample flagging.
- SESAME_POOBAH_QC: SeSAMe pOOBAH-based detection p-value filtering; outputs normalised beta matrix and filtered RGset.
- QC_REPORT_CONTROLS: minfi control probe QC report and bisulphite conversion efficiency plots.
- SNP_HEATMAP: SNP-based sample identity and relatedness heatmap.
- DENSITY_PLOTS: Beta-value density plots for raw and SeSAMe-normalised data.
- SEX_QC: Predicted sex computation and sex-discordance flagging.
- FILTER_XY_NONCG_BEADS: Independent filtering of sex-chromosome probes and non-CpG probes.
- CELL_COMP: Optional blood cell composition estimation using FlowSorted.Blood.EPIC (enabled with `--do_estimate_cellcomp`).
- SPLIT_COLLAPSE: SNP/cross-hybridisation probe removal and EPICv2 positional-replicate collapsing; produces separate beta matrices for DMP and DMR analyses.
- FINAL_EXPORTS: Export of final beta and M-value matrices as RDS (and optionally CSV).
- ALIGN_META: Alignment of methylation matrices with clinical metadata; optional ComBat batch correction.
- COVARIATE_PCA: Principal component analysis of beta values versus clinical covariates; optional ChAMP SVD analysis (enabled with `--run_covariate_pca`).
- DMP_LIMMA: Differential methylation position analysis using limma with optional ComBat batch correction (enabled with `--run_dmp`).
- DMR_DMRCATE: Differential methylation region analysis using DMRcate with EPICv2 replicate remapping (enabled with `--run_dmr`).
- LOG_VERSIONS: Logging of all R/Bioconductor package versions for reproducibility.
- Parameters: `--annotation`, `--genome_build`, `--p_det`, `--fail_sample_fr`, `--fail_probe_fr`, `--sesame_prep`, `--do_rmsnpandch`, `--snp_dist_bp`, `--snp_maf`, `--do_estimate_cellcomp`, `--filter_xy`, `--filter_noncg`, `--export_csv_matrices`, `--save_final_matrices_rds`, `--matrix_rds_compression`, `--keep_groups`, `--epicv2_remap`, `--posreps_strategy`, `--run_covariate_pca`, `--run_champ_svd`, `--run_dmp`, `--run_dmr`.
- Bead count masking (minimum 3 beads per probe) is applied internally by SeSAMe's `qualityMask` (Q step) inside `SESAME_POOBAH_QC`, replacing the former two-stage approach that also applied `minfi::getNBeads()` in `FILTER_XY_NONCG_BEADS`. The parameters `--filter_beads`, `--bead_min`, and `--bead_fail_frac` have been removed.
- Support for Docker, Singularity, Conda/Mamba, Podman, Charliecloud, Apptainer and Wave execution profiles.
- nf-core template v3.5.2.

### `Fixed`

### `Dependencies`

- R 4.4.0 (Bioconductor 3.19 base image)
- minfi 1.48.0
- sesame 1.18.0
- sesameData 1.18.0
- limma 3.58.0
- DMRcate 2.16.0
- FlowSorted.Blood.EPIC 2.4.0
- ChAMP 2.32.0
- sva 3.50.0
- missMethyl 1.36.0
- pheatmap 1.0.12
- ggplot2 3.5.1
- dplyr 1.1.4
- readr 2.1.5
- optparse 1.7.5
- yaml 2.3.10

### `Deprecated`

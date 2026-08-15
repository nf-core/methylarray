# nf-core/methylarray: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/methylarray/usage](https://nf-co.re/methylarray/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/methylarray processes Illumina EPICv2 DNA methylation array data starting from raw IDAT files. It requires three inputs: a sample sheet (TSV), a directory of IDAT files, and a metadata CSV with sample phenotype or clinical information. This page describes how to prepare those inputs and how to run the pipeline.

## Samplesheet input (`--input`)

The sample sheet must be a **tab-separated** file (`.tsv`) with a header row and exactly three columns. Each row describes one sample.

```bash
--input '[path to samplesheet.tsv]'
```

```tsv title="samplesheet.tsv"
Sample_Name	Sentrix_ID	Sentrix_Position
SAMPLE_001	207842290093	R01C01
SAMPLE_002	207842290093	R01C02
SAMPLE_003	207842290094	R01C01
```

| Column             | Description                                                                                                        |
| ------------------ | ------------------------------------------------------------------------------------------------------------------ |
| `Sample_Name`      | Unique sample identifier. Must match the `sample_id` column in the metadata CSV.                                   |
| `Sentrix_ID`       | Illumina array barcode (e.g. `207842290093`). Together with `Sentrix_Position` this identifies the IDAT file pair. |
| `Sentrix_Position` | Array section identifier in the format `RxxCxx` (e.g. `R01C01`).                                                   |

The pipeline locates IDAT files by concatenating `Sentrix_ID` and `Sentrix_Position` (e.g. `207842290093_R01C01_Grn.idat` and `207842290093_R01C01_Red.idat`). Files may reside in subdirectories within `--idat_dir`.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## IDAT directory (`--idat_dir`)

Point this parameter to the directory containing your raw IDAT files. Subdirectories are searched recursively.

```bash
--idat_dir '[path to directory containing IDAT files]'
```

Each sample requires a green-channel (`*_Grn.idat`) and a red-channel (`*_Red.idat`) file. Files are matched to samplesheet rows using the `Sentrix_ID_Sentrix_Position` prefix.

## Metadata file (`--meta_file`)

A **comma-separated** (CSV) or tab-separated file providing sample-level phenotype and clinical information. The file must include a header row. The `sample_id` column must match the `Sample_Name` values in the samplesheet — only samples present in both files will be carried forward into differential methylation analysis.

```bash
--meta_file '[path to metadata.csv]'
```

### Required columns

| Column      | Type   | Description                                                                                                                                               |
| ----------- | ------ | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample_id` | string | Must match `Sample_Name` in the samplesheet exactly.                                                                                                      |
| `group`     | string | Categorical group variable used as the primary factor in DMP/DMR models. The level `"CONTROL"` is used as the reference level for all pairwise contrasts. |
| `sex`       | string | Reported biological sex (e.g. `"M"` / `"F"`). Used to validate predicted sex from methylation data.                                                       |

### Optional columns

| Column  | Type    | Description                                                                                                                                                         |
| ------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `age`   | numeric | Sample age. If present, included as a covariate in DMP/DMR linear models.                                                                                           |
| `bmi`   | numeric | Body mass index. If present, included as a covariate in DMP/DMR linear models.                                                                                      |
| `Plate` | string  | Batch identifier (e.g. hybridization plate). Used for ComBat batch correction when the `DMP_limma__combat` and `DMR_DMRcate_EPICv2__combat` analysis modes are run. |
| `CD8T`  | numeric | Estimated CD8+ T cell proportion (from cell composition step or user-supplied).                                                                                     |
| `CD4T`  | numeric | Estimated CD4+ T cell proportion.                                                                                                                                   |
| `NK`    | numeric | Estimated NK cell proportion.                                                                                                                                       |
| `Bcell` | numeric | Estimated B cell proportion.                                                                                                                                        |
| `Mono`  | numeric | Estimated monocyte proportion.                                                                                                                                      |
| `Neu`   | numeric | Estimated neutrophil proportion.                                                                                                                                    |
| `Gran`  | numeric | Estimated granulocyte proportion.                                                                                                                                   |

Cell composition columns are appended automatically when `--do_estimate_cellcomp true` (the default). If you pre-estimated cell proportions, you can supply them directly in the metadata CSV and disable estimation with `--do_estimate_cellcomp false`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/methylarray \
   --input ./samplesheet.tsv \
   --idat_dir ./idat_files/ \
   --meta_file ./metadata.csv \
   --outdir ./results \
   -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/methylarray -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/methylarray
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/methylarray releases page](https://github.com/nf-core/methylarray/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

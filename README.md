<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-ontvar_logo_dark.png">
    <img alt="nf-core/ontvar" src="docs/images/nf-core-ontvar_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/ontvar/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/ontvar/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/ontvar/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/ontvar/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/ontvar/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/ontvar)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23ontvar-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/ontvar)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/ontvar** is a comprehensive structural variant calling, filtering, annotation and consensus set generation pipeline for ONT's long read sequencing data.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Workflow Overview

```mermaid 
---
config:
  layout: elk
---
flowchart TB
    INPUT_SAMPLE["Sample Sheet"] L_INPUT_SAMPLE_SPLIT_0@--> SPLIT["Split by Type"]
    INPUT_ANNOT["AnnotSV Annotations"] L_INPUT_ANNOT_ANNOT_SETUP_0@--> ANNOT_SETUP["Annotations<br>Available?"]
    SPLIT L_SPLIT_CASES_0@--> CASES["Case Samples"] & CONTROLS["Control Samples"]
    CASES L_CASES_SNIFFLES_0@--> SNIFFLES["Sniffles"] & CUTESV["cuteSV"] & SEVERUS_CTRL["Severus<br>(with control)"] & SEVERUS_NO["Severus<br>(no control)"]
    CONTROLS L_CONTROLS_SEVERUS_CTRL_0@--> SEVERUS_CTRL
    SNIFFLES L_SNIFFLES_FMT1_0@--> FMT1["Format VCF Headers"]
    CUTESV L_CUTESV_FMT2_0@--> FMT2["Format VCF Headers"]
    SEVERUS_CTRL L_SEVERUS_CTRL_FMT3_0@--> FMT3["Format VCF Headers"]
    SEVERUS_NO L_SEVERUS_NO_FMT3_0@--> FMT3
    FMT1 L_FMT1_RAW_VCFS_0@--> RAW_VCFS["Individual Caller VCFs"]
    FMT2 L_FMT2_RAW_VCFS_0@--> RAW_VCFS
    FMT3 L_FMT3_RAW_VCFS_0@--> RAW_VCFS
    RAW_VCFS L_RAW_VCFS_SUM_RAW_0@--> SUM_RAW["Summarize Raw Calls"] & JASMINE_SAMPLE["Jasmine: Merge Callers<br>(per sample)"]
    SUM_RAW L_SUM_RAW_PLOT_RAW_0@--> PLOT_RAW["Plot: Raw Caller Counts"]
    JASMINE_SAMPLE L_JASMINE_SAMPLE_FIX_HDR_0@--> FIX_HDR["Fix VCF Headers"]
    FIX_HDR L_FIX_HDR_FILTER_CHR_0@--> FILTER_CHR["Filter Chromosomes"]
    FILTER_CHR L_FILTER_CHR_SORT1_0@--> SORT1["Sort VCF"]
    SORT1 --> SUPPORT_FILTER["Filter: ≥2 Caller Support"]
    SUPPORT_FILTER L_SUPPORT_FILTER_SUM_CONSENSUS_0@--> SUM_CONSENSUS["Summarize Consensus Calls"] & SVDB_SAMPLE["Annotate: Population AF<br>(sample)"] & COLLECT["Sample VCFs"]
    SUM_CONSENSUS L_SUM_CONSENSUS_PLOT_CONS_0@--> PLOT_CONS["Plot: Consensus Counts"]
    SVDB_SAMPLE L_SVDB_SAMPLE_ANNOT_SAMPLE_RAW_0@--> ANNOT_SAMPLE_RAW["AnnotSV: FUNCTIONAL ANNOTATIONS<br>(sample: before AF filter)"] & AF_SAMPLE["Filter: Allele Frequency<br>(per sample)"]
    AF_SAMPLE L_AF_SAMPLE_SUM_FILT_0@--> SUM_FILT["Summarize Filtered Calls"] & ANNOT_SAMPLE_FINAL["AnnotSV: FUNCTIONAL ANNOTATIONS <br>(sample: after AF FILTER)"]
    SUM_FILT L_SUM_FILT_PLOT_FILT_0@--> PLOT_FILT["Plot: Filtered Counts"]
    COLLECT L_COLLECT_JASMINE_COHORT_0@--> JASMINE_COHORT["Jasmine: Merge Samples<br>(cohort-wide)"]
    JASMINE_COHORT L_JASMINE_COHORT_SVDB_COHORT_0@--> SVDB_COHORT["Annotate: Population AF<br>(cohort)"]
    SVDB_COHORT L_SVDB_COHORT_ANNOT_COHORT_RAW_0@--> ANNOT_COHORT_RAW["AnnotSV: FUNCTIONAL ANNOTATIONS<br>(cohort: before AF filter)"] & SUM_COHORT_RAW["Summarize Cohort<br>(annotated)"] & AF_COHORT["Filter: Allele Frequency<br>(cohort)"]
    SUM_COHORT_RAW L_SUM_COHORT_RAW_PLOT_COHORT_RAW_0@--> PLOT_COHORT_RAW["Plot: Cohort Annotated"]
    AF_COHORT L_AF_COHORT_SUM_COHORT_FILT_0@--> SUM_COHORT_FILT["Summarize Cohort<br>(filtered)"] & ANNOT_COHORT_FINAL["AnnotSV: FUNCTIONAL ANNOTATIONS<br>(cohort: after AF filter)"]
    SUM_COHORT_FILT L_SUM_COHORT_FILT_PLOT_COHORT_FILT_0@--> PLOT_COHORT_FILT["Plot: Cohort Filtered"]
    ANNOT_SETUP L_ANNOT_SETUP_INSTALL_0@-- No --> INSTALL["Install AnnotSV data"]
    ANNOT_SETUP L_ANNOT_SETUP_CHECK_TAR_0@-- Yes --> CHECK_TAR["tar.gz?"]
    CHECK_TAR L_CHECK_TAR_UNTAR_0@-- Yes --> UNTAR["Extract Archive"]
    CHECK_TAR L_CHECK_TAR_ANNOT_DB_0@-- No --> ANNOT_DB["AnnotSV data"]
    INSTALL L_INSTALL_ANNOT_DB_0@--> ANNOT_DB
    UNTAR L_UNTAR_ANNOT_DB_0@--> ANNOT_DB
    ANNOT_DB L_ANNOT_DB_ANNOT_SAMPLE_RAW_0@--> ANNOT_SAMPLE_RAW & ANNOT_SAMPLE_FINAL & ANNOT_COHORT_RAW & ANNOT_COHORT_FINAL
    ANNOT_SAMPLE_FINAL L_ANNOT_SAMPLE_FINAL_OUT_SAMPLE_0@--> OUT_SAMPLE["Sample-Level Results"]
    ANNOT_COHORT_FINAL L_ANNOT_COHORT_FINAL_OUT_COHORT_0@--> OUT_COHORT["Cohort-Level Results"]
    INPUT_REF["Reference Genome"]
    INPUT_SAMPLE@{ shape: rounded}
    SPLIT@{ shape: rounded}
    INPUT_ANNOT@{ shape: rounded}
    ANNOT_SETUP@{ shape: rounded}
    CASES@{ shape: rounded}
    CONTROLS@{ shape: rounded}
    SNIFFLES@{ shape: rounded}
    CUTESV@{ shape: rounded}
    SEVERUS_CTRL@{ shape: rounded}
    SEVERUS_NO@{ shape: rounded}
    FMT1@{ shape: rounded}
    FMT2@{ shape: rounded}
    FMT3@{ shape: rounded}
    RAW_VCFS@{ shape: rounded}
    SUM_RAW@{ shape: rounded}
    JASMINE_SAMPLE@{ shape: rounded}
    PLOT_RAW@{ shape: rounded}
    FIX_HDR@{ shape: rounded}
    FILTER_CHR@{ shape: rounded}
    SORT1@{ shape: rounded}
    SUPPORT_FILTER@{ shape: rounded}
    SUM_CONSENSUS@{ shape: rounded}
    SVDB_SAMPLE@{ shape: rounded}
    COLLECT@{ shape: rounded}
    PLOT_CONS@{ shape: rounded}
    ANNOT_SAMPLE_RAW@{ shape: rounded}
    AF_SAMPLE@{ shape: rounded}
    SUM_FILT@{ shape: rounded}
    ANNOT_SAMPLE_FINAL@{ shape: rounded}
    PLOT_FILT@{ shape: rounded}
    JASMINE_COHORT@{ shape: rounded}
    SVDB_COHORT@{ shape: rounded}
    ANNOT_COHORT_RAW@{ shape: rounded}
    SUM_COHORT_RAW@{ shape: rounded}
    AF_COHORT@{ shape: rounded}
    PLOT_COHORT_RAW@{ shape: rounded}
    SUM_COHORT_FILT@{ shape: rounded}
    ANNOT_COHORT_FINAL@{ shape: rounded}
    PLOT_COHORT_FILT@{ shape: rounded}
    INSTALL@{ shape: rounded}
    CHECK_TAR@{ shape: rounded}
    UNTAR@{ shape: rounded}
    ANNOT_DB@{ shape: rounded}
    OUT_SAMPLE@{ shape: rounded}
    OUT_COHORT@{ shape: rounded}
    INPUT_REF@{ shape: rounded}
     INPUT_SAMPLE:::input
     SPLIT:::decision
     INPUT_ANNOT:::input
     ANNOT_SETUP:::decision
     SNIFFLES:::caller
     CUTESV:::caller
     SEVERUS_CTRL:::caller
     SEVERUS_NO:::caller
     FMT1:::process
     FMT2:::process
     FMT3:::process
     SUM_RAW:::summary
     JASMINE_SAMPLE:::process
     PLOT_RAW:::output
     FIX_HDR:::process
     FILTER_CHR:::process
     SORT1:::process
     SUPPORT_FILTER:::filter
     SUM_CONSENSUS:::summary
     SVDB_SAMPLE:::process
     PLOT_CONS:::output
     ANNOT_SAMPLE_RAW:::process
     AF_SAMPLE:::filter
     SUM_FILT:::summary
     ANNOT_SAMPLE_FINAL:::process
     PLOT_FILT:::output
     JASMINE_COHORT:::process
     SVDB_COHORT:::process
     ANNOT_COHORT_RAW:::process
     SUM_COHORT_RAW:::summary
     AF_COHORT:::filter
     PLOT_COHORT_RAW:::output
     SUM_COHORT_FILT:::summary
     ANNOT_COHORT_FINAL:::process
     PLOT_COHORT_FILT:::output
     INSTALL:::process
     CHECK_TAR:::decision
     UNTAR:::process
     OUT_SAMPLE:::output
     OUT_COHORT:::output
     INPUT_REF:::input
    classDef input fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    classDef caller fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    classDef process fill:#e8f5e9,stroke:#1b5e20,stroke-width:2px
    classDef filter fill:#fff3e0,stroke:#e65100,stroke-width:2px
    classDef summary fill:#fff8e1,stroke:#f57f17,stroke-width:2px
    classDef output fill:#c8e6c9,stroke:#2e7d32,stroke-width:3px
    classDef decision fill:#fce4ec,stroke:#880e4f,stroke-width:2px
    L_INPUT_SAMPLE_SPLIT_0@{ animation: fast } 
    L_INPUT_ANNOT_ANNOT_SETUP_0@{ animation: fast } 
    L_SPLIT_CASES_0@{ animation: fast } 
    L_SPLIT_CONTROLS_0@{ animation: fast } 
    L_CASES_SNIFFLES_0@{ animation: fast } 
    L_CASES_CUTESV_0@{ animation: fast } 
    L_CASES_SEVERUS_NO_0@{ animation: fast } 
    L_CONTROLS_SEVERUS_CTRL_0@{ animation: fast } 
    L_SNIFFLES_FMT1_0@{ animation: fast } 
    L_CUTESV_FMT2_0@{ animation: fast } 
    L_SEVERUS_CTRL_FMT3_0@{ animation: fast } 
    L_SEVERUS_NO_FMT3_0@{ animation: fast } 
    L_FMT1_RAW_VCFS_0@{ animation: fast } 
    L_FMT2_RAW_VCFS_0@{ animation: fast } 
    L_FMT3_RAW_VCFS_0@{ animation: fast } 
    L_RAW_VCFS_SUM_RAW_0@{ animation: fast } 
    L_RAW_VCFS_JASMINE_SAMPLE_0@{ animation: fast } 
    L_SUM_RAW_PLOT_RAW_0@{ animation: fast } 
    L_JASMINE_SAMPLE_FIX_HDR_0@{ animation: fast } 
    L_FIX_HDR_FILTER_CHR_0@{ animation: fast } 
    L_FILTER_CHR_SORT1_0@{ animation: fast } 
    L_SUPPORT_FILTER_SUM_CONSENSUS_0@{ animation: none } 
    L_SUPPORT_FILTER_SVDB_SAMPLE_0@{ animation: fast } 
    L_SUPPORT_FILTER_COLLECT_0@{ animation: fast } 
    L_SUM_CONSENSUS_PLOT_CONS_0@{ animation: none } 
    L_SVDB_SAMPLE_ANNOT_SAMPLE_RAW_0@{ animation: fast } 
    L_SVDB_SAMPLE_AF_SAMPLE_0@{ animation: fast } 
    L_AF_SAMPLE_SUM_FILT_0@{ animation: fast } 
    L_AF_SAMPLE_ANNOT_SAMPLE_FINAL_0@{ animation: fast } 
    L_SUM_FILT_PLOT_FILT_0@{ animation: fast } 
    L_COLLECT_JASMINE_COHORT_0@{ animation: fast } 
    L_JASMINE_COHORT_SVDB_COHORT_0@{ animation: fast } 
    L_SVDB_COHORT_ANNOT_COHORT_RAW_0@{ animation: fast } 
    L_SVDB_COHORT_SUM_COHORT_RAW_0@{ animation: fast } 
    L_SVDB_COHORT_AF_COHORT_0@{ animation: fast } 
    L_SUM_COHORT_RAW_PLOT_COHORT_RAW_0@{ animation: fast } 
    L_AF_COHORT_SUM_COHORT_FILT_0@{ animation: fast } 
    L_AF_COHORT_ANNOT_COHORT_FINAL_0@{ animation: fast } 
    L_SUM_COHORT_FILT_PLOT_COHORT_FILT_0@{ animation: fast } 
    L_ANNOT_SETUP_INSTALL_0@{ animation: slow } 
    L_ANNOT_SETUP_CHECK_TAR_0@{ animation: slow } 
    L_CHECK_TAR_UNTAR_0@{ animation: fast } 
    L_CHECK_TAR_ANNOT_DB_0@{ animation: fast } 
    L_INSTALL_ANNOT_DB_0@{ animation: fast } 
    L_UNTAR_ANNOT_DB_0@{ animation: fast } 
    L_ANNOT_DB_ANNOT_SAMPLE_RAW_0@{ animation: fast } 
    L_ANNOT_DB_ANNOT_SAMPLE_FINAL_0@{ animation: fast } 
    L_ANNOT_DB_ANNOT_COHORT_RAW_0@{ animation: fast } 
    L_ANNOT_DB_ANNOT_COHORT_FINAL_0@{ animation: fast } 
    L_ANNOT_SAMPLE_FINAL_OUT_SAMPLE_0@{ animation: fast } 
    L_ANNOT_COHORT_FINAL_OUT_COHORT_0@{ animation: fast }

```

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/ontvar \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/ontvar/usage) and the [parameter documentation](https://nf-co.re/ontvar/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/ontvar/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/ontvar/output).

## Credits

nf-core/ontvar was originally written by Manas Sehgal.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#ontvar` channel](https://nfcore.slack.com/channels/ontvar) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/ontvar for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

# MACE Metagenomics: Snakemake workflows for metagenomic assembly, classification, binning, and MAG summarisation

This repository contains a modular Snakemake-based pipeline for metagenomic assembly, contig classification, genome binning, and downstream summarisation of both metagenomes and MAGs.

The pipeline is designed around soil metagenomes and supports:
- **co-assembly** of all reads together
- **sub-co-assembly** of groups of samples
- **single-sample assembly**

The main workflow files in this branch are:
- `metagenome_assemble.smk`
- `metagenome_classify.smk`
- `metagenome_binning.smk`
- `summarise_metagenome.smk`
- `summarise_mags.smk`

---

## Overview

### 1. Assembly: `metagenome_assemble.smk`
This workflow performs metagenomic assembly starting from short reads and, optionally, long reads for hybrid assemblies.

In brief:
1. Short reads are trimmed with **trim_galore**
2. If long reads are provided, they can be preprocessed with **porechop** and/or **chopper**
3. Reads are assembled into contigs using:
   - **MEGAHIT** for short-read assemblies
   - **MetaSPAdes** for hybrid assemblies

The short reads are expected to be FASTQ files in the input short-read folder, using sample naming compatible with:
- `*_R1_001.f*q.gz`
- `*_R2_001.f*q.gz`

Assembly can be run in three modes:
- `co`: all reads assembled together into one co-assembly
- `sub`: reads assembled into several grouped sub-co-assemblies based on the metadata `Group` column
- `single`: each sample assembled independently

### 2. Contig classification: `metagenome_classify.smk`
This workflow separates contigs into **prokaryotic**, **eukaryotic**, and **other** categories.

It is used to focus downstream binning on the **prokaryotic contigs**.

At present this workflow uses:
- **WhoKaryote**
- **Tiara**

The key downstream output is:
- `prokaryotic_contigs.fasta`

### 3. Binning: `metagenome_binning.smk`
This workflow performs MAG recovery from the prokaryotic contigs.

The current approach is:
1. Filter prokaryotic contigs
2. Rename contigs
3. Compute coverage information with **CoverM**
4. Run three binners:
   - **MetaBAT2**
   - **ComeBin**
   - **MaxBin2**
5. Combine and optimise bins with **Binette**

The result is a set of bins/MAGs for each assembly suffix.

---

## Downstream summarisation workflows

### 4. Metagenome-level summary: `summarise_metagenome.smk`
This workflow summarises the metagenome assemblies themselves.

It currently produces:
- **CoverM** output for **trimmed_mean** abundance of contigs across samples
- functional annotation using **eggNOG-mapper**

This is done per assembly.

### 5. MAG-level summary: `summarise_mags.smk`
This workflow summarises the recovered MAGs.

It currently performs:
- **CheckM2** to assess completeness and contamination
- **Galah** to dereplicate MAGs
- **Bakta** for functional annotation
- **CoverM** for MAG abundance across all samples
- **GTDB-Tk** for taxonomy assignment

The MAG abundance step uses the trimmed reads found under the output directory, and the MAG workflow uses all assemblies listed in the `suffixes` config entry.

---

## Important pipeline defaults

### MEGAHIT preset
For short-read assemblies, **MEGAHIT** currently uses the **`meta-large`** preset because this pipeline is mainly used on **soil metagenomes**.

This is currently hard-coded in the Snakefile. In the future it would be better to expose this as a config parameter.

### Galah dereplication threshold
For MAG dereplication, **Galah** is currently run at **99% ANI**, i.e. approximately **strain-level dereplication**.

This is also currently hard-coded and would ideally become a user-settable parameter in the future.

---

## Configuration

At present the repository contains separate config files:
- `config_assemble.yml`
- `config_binning.yml`

A merged config is possible and likely cleaner long term, but for now the workflows still reflect separate assembly and binning configuration.

### Current assembly config example
```yaml
name: my_dataset
illumina_folder: /data/metagenomics/raw_data/path/to/folder/with/fastqfiles
nanopore_folder: null   # set to a folder to enable hybrid assemblies
metadata: path/to/final_metadata.tsv

assembly_type: "sub"    # can also include "co" and/or "single"
cpus: 64

# preprocessing
trim_galore_threads: 8
long_reads_preprocessing:
  porechop: true
  chopper: true

# output
outdir: /data/metagenomics/processed_data/name_of_the_output_results
```

### Current binning config example
```yaml
outdir: "/data/metagenomics/processed_data/my_results"
suffixes: ["sub_assembly_1", "sub_assembly_2", "sub_assembly_3"]
threads: 32
short_reads_folder: "/data/metagenomics/processed_data/my_results/trimmed_reads/short_reads"
```

### Suggested merged config example
A merged config would make the relationship between assembly, binning, and summarisation clearer. The metadata.tsv file is used for sub-co-assemblies to specify in the column `Group` of a .tsv file which samples belong to it in the column `Sample`.

```yaml
name: my_dataset

# inputs
illumina_folder: /data/metagenomics/raw_data/path/to/folder/with/fastqfiles
nanopore_folder: null
metadata: path/to/metadata.tsv

# assembly mode
assembly_type: ["sub"]   # choose from "co", "sub", "single"

# compute
threads: 64
trim_galore_threads: 8

# long-read preprocessing
long_reads_preprocessing:
  porechop: true
  chopper: true

# output
outdir: /data/metagenomics/processed_data/name_of_the_output_results

# downstream binning / MAG summary
suffixes: [sub_assembly_1, sub_assembly_2, sub_assembly_3]
```

### Notes on config cleanup
A few config-related improvements are still needed:
- the repository currently mixes `cpus` and `threads`; these should ideally be standardised
- if configs are merged, all workflow references should consistently use one naming convention, preferably `threads`
- `short_reads_folder` should ideally be derived automatically from `outdir` as:
  - `{outdir}/trimmed_reads/short_reads`

---

## Running the workflows

### 1. Assembly
```bash
snakemake --snakefile metagenome_assemble.smk --configfile config_assemble.yml --use-conda --cores 32
```

### 2. Classification
This workflow expects the assembly outputs, especially `assembly_groups.tsv`, to already exist in the output directory.

```bash
snakemake --snakefile metagenome_classify.smk --config outdir=/path/to/results cpus=8 --use-conda --cores 32
```

### 3. Binning
```bash
snakemake --snakefile metagenome_binning.smk --configfile config_binning.yml --use-conda --cores 32
```

### 4. Metagenome summarisation
```bash
snakemake --snakefile summarise_metagenome.smk --configfile config_binning.yml --use-conda --cores 32
```

### 5. MAG summarisation
```bash
snakemake --snakefile summarise_mags.smk --configfile config_binning.yml --use-conda --cores 32
```

---

## Conda environments

Environment YAMLs in this repository are stored in:
- `envs/`

Examples present in this branch include:
- `envs/MACE_metagenomics.yml`
- `envs/binette.yml`
- `envs/binning_base.yml`
- `envs/checkm2.yml`
- `envs/comebin.yml`
- `envs/coverm.yml`
- `envs/gtdbtk-2.7.2.yml`
- `envs/maxbin2.yml`
- etc.

So at the moment, some “envs” are actually referenced as **environment names** rather than **YAML files**.

## Expected folder structure

```text
MACE_metagenomics/
├── README.md
├── config_assemble.yml
├── config_binning.yml
├── metagenome_assemble.smk
├── metagenome_classify.smk
├── metagenome_binning.smk
├── summarise_metagenome.smk
├── summarise_mags.smk
└── envs/
    ├── MACE_metagenomics.yml
    ├── binette.yml
    ├── binning_base.yml
    ├── checkm2.yml
    ├── comebin.yml
    ├── coverm.yml
    ├── gtdbtk-2.7.2.yml
    ├── maxbin2.yml
    ├── metabat2.yml
    ├── semibin2.yml
    ├── whokaryote.yml
    └── ...
```

---

## Notes and caveats

- The classification workflow currently depends on `assembly_groups.tsv` already existing in the output directory from the assembly step.
- Some environment references should be cleaned up so all rules consistently point to `envs/*.yml`.
- Some parameters that are currently hard-coded in the Snakefiles would be better exposed in config files:
  - MEGAHIT preset
  - Galah ANI threshold
  - database paths
  - some thread / CPU settings
- `short_reads_folder` should ideally not need to be manually repeated if it can be derived from `outdir`.

---

## Troubleshooting

- Always use `--use-conda`
- Check that your sample names in the metadata match the FASTQ names
- Check that reads follow the expected naming convention:
  - `sample_R1_001.fastq.gz`
  - `sample_R2_001.fastq.gz`
- Confirm that `metadata` contains the required columns, especially `Sample`, and for sub-assemblies also `Group`
- If classification or summarisation fails, verify that upstream outputs already exist in `outdir`
- If a rule fails because of a missing environment, check whether the workflow expects:
  - an explicit `.yml` file in `envs/`, or
  - a pre-existing named conda environment

---

## Getting help

For bugs, missing environment files, or workflow improvements, please open an issue in this repository.

This pipeline is still being cleaned up and documented, so issues and suggested improvements are welcome.

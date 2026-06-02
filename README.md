# MACE metagenomics pipeline

This repository contains Snakemake workflows for metagenomic assembly, contig classification, binning, and metagenome summarisation.

## Workflows

- `metagenome_assemble.smk`: trims reads with Trim Galore and assembles contigs with MEGAHIT, or MetaSPAdes for hybrid assemblies.
- `metagenome_classify.smk`: separates eukaryotic/prokaryotic/other contigs and carries prokaryotic contigs forward.
- `metagenome_binning.smk`: bins prokaryotic contigs using MetaBAT2, COMEbin, and MaxBin2, then combines/optimizes bins with Binette.
- `summarise_metagenome.smk`: produces contig abundance summaries using CoverM `trimmed_mean` across samples for each assembly and functional annotation using eggNOG-mapper.
- `summarise_mags.smk`: MAG-level summary workflow (CheckM2 completeness/contamination, Galah dereplication, Bakta annotation, MAG abundance across all samples from assemblies listed in `suffixes`, and GTDB-Tk taxonomy).

## Assembly modes

Set in `assembly.assembly_type`:

- `co`: all reads together
- `sub`: subgroup co-assemblies based on metadata
- `single`: per-sample assemblies

## Hardcoded settings

- MEGAHIT currently uses the `meta-large` preset (soil-focused default in this pipeline).
- Galah dereplication is currently hardcoded at 99% ANI (strain-level).

These are currently hardcoded and can be made configurable in future updates.

## Config

The workflows use one merged config file (`test_set/config.yaml`) with shared sections (`inputs`, `resources`, `assembly`, `binning`).

`binning.short_reads_folder` is no longer required in config; workflows derive it from:

- `{outdir}/trimmed_reads/short_reads`

### Example merged config

```yaml
name: test_set
outdir: /data/metagenomics/processed_data/test_set

inputs:
  illumina_folder: /data/metagenomics/raw_data/test_data
  nanopore_folder: null
  metadata: /data/metagenomics/raw_data/test_data/test_metadata.csv

resources:
  threads: 64
  trim_galore_threads: 8

assembly:
  assembly_type:
    - sub
  long_reads_preprocessing:
    porechop: true
    chopper: true

binning:
  suffixes:
    - sub_assembly_Group1
```

## Run order

```bash
snakemake -s metagenome_assemble.smk --use-conda --configfile test_set/config.yaml -j 64
snakemake -s metagenome_classify.smk --use-conda --configfile test_set/config.yaml -j 64
snakemake -s metagenome_binning.smk --use-conda --configfile test_set/config.yaml -j 64
snakemake -s summarise_metagenome.smk --use-conda --configfile test_set/config.yaml -j 64
```

## Conda environments

Workflows use YAML-backed envs under `envs/`, including:

- `envs/comebin.yml`
- `envs/checkm2.yml`
- `envs/gtdbtk-2.7.2.yml`
- `envs/bakta.yml`
- `envs/eggnog_blast.yml`

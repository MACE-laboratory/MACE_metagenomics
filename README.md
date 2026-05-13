# MACE Metagenomics: Snakemake Workflows

Welcome to the MACE Metagenomics repository.
This project provides modular, scalable, and reproducible Snakemake workflows for metagenomic assembly, classification, binning, and downstream summarisation.

The workflow has now been simplified to use **one shared config file per dataset/run**.
For example, the test dataset now uses:

`test_set/config.yaml`

---

## 📂 Available Workflows

| Workflow | Purpose |
|---|---|
| `metagenome_assemble_march2026.smk` | Metagenome assembly from short and optional long reads |
| `metagenome_classify_march2026.smk` | Contig classification into prokaryotic / eukaryotic / other |
| `metagenome_binning_march2026.smk` | Genome binning from classified prokaryotic contigs |
| `summarise_metagenome_v01.smk` | Functional and coverage summarisation of classified assemblies |
| `summarise_bins.smk` | Coming soon: dereplication, taxonomy, quality, and metabolic profiling of bins |

All workflows are intended to run with the same unified config file.

---

## ⚙️ One-config workflow

The pipeline now uses a single YAML config with shared settings plus step-specific sections.

### Example unified config

````yaml name=test_set/config.yaml url=https://github.com/MACE-laboratory/MACE_metagenomics/blob/v01_dev/test_set/config.yaml
name: test_set
outdir: /data/metagenomics/processed_data/test_set

inputs:
  illumina_folder: /data/metagenomics/raw_data/test_data
  nanopore_folder: null
  metadata: /data/metagenomics/raw_data/test_data/test_metadata.csv

resources:
  cpus: 64
  threads: 64
  trim_galore_threads: 8

assembly:
  assembly_type:
    - sub
  long_reads_preprocessing:
    porechop: true
    chopper: true

binning:
  short_reads_folder: /data/metagenomics/processed_data/test_set/trimmed_reads/short_reads
  suffixes:
    - sub_assembly_Group1
````

### Config structure

| Section | Purpose |
|---|---|
| `name` | Name of the run or dataset |
| `outdir` | Base output directory for all steps |
| `inputs` | Input data locations such as Illumina, Nanopore, and metadata |
| `resources` | Shared CPU/thread settings |
| `assembly` | Assembly-specific parameters |
| `binning` | Binning-specific parameters |

### Notes

- `inputs.nanopore_folder: null` means short-read-only assembly.
- `assembly.assembly_type` can contain one or more of:
  - `co`
  - `sub`
  - `single`
- `binning.suffixes` can be set explicitly, or the workflow can infer assemblies from `assembly_groups.tsv`.

---

## 🔬 Workflow tutorial: run the full pipeline

This section describes the recommended order for running the whole metagenomics workflow:

1. assembly
2. classify
3. binning
4. summarise
5. summarise bins (coming soon)

In all examples below, the same config file is used:

```bash name=commands_config.sh
CONFIG=test_set/config.yaml
```

### 1. Assembly

**Workflow:** `metagenome_assemble_march2026.smk`

This step:
- validates metadata against the raw read files
- trims Illumina reads
- optionally preprocesses Nanopore reads
- creates co-assembly, subgroup assembly, and/or single-sample assembly outputs depending on config
- writes `assembly_groups.tsv`

#### Dry run
```bash name=assembly_dry_run.sh
snakemake -s metagenome_assemble_march2026.smk --use-conda --configfile test_set/config.yaml -n
```

#### Run
```bash name=assembly_run.sh
snakemake -s metagenome_assemble_march2026.smk --use-conda --configfile test_set/config.yaml -j 64
```

#### Main outputs
- `{outdir}/assemblies/`
- `{outdir}/assembly_groups.tsv`
- `{outdir}/trimmed_reads/`
- `{outdir}/assemblies/tools_versions.yaml`

---

### 2. Contig classification

**Workflow:** `metagenome_classify_march2026.smk`

This step:
- uses assembly outputs from the previous step
- classifies contigs with Tiara / WhoKaryote-related logic in the workflow
- splits assemblies into:
  - `prokaryotic_contigs.fasta`
  - `eukaryotic_contigs.fasta`
  - `other_contigs.fasta`

#### Dry run
```bash name=classify_dry_run.sh
snakemake -s metagenome_classify_march2026.smk --use-conda --configfile test_set/config.yaml -n
```

#### Run
```bash name=classify_run.sh
snakemake -s metagenome_classify_march2026.smk --use-conda --configfile test_set/config.yaml -j 64
```

#### Important
This workflow expects `assembly_groups.tsv` and assembly outputs to already exist, so run assembly first.

---

### 3. Genome binning

**Workflow:** `metagenome_binning_march2026.smk`

This step:
- filters and renames prokaryotic contigs
- computes coverage with CoverM
- runs multiple binning tools
- consolidates bins with Binette

#### Dry run
```bash name=binning_dry_run.sh
snakemake -s metagenome_binning_march2026.smk --use-conda --configfile test_set/config.yaml -n
```

#### Run
```bash name=binning_run.sh
snakemake -s metagenome_binning_march2026.smk --use-conda --configfile test_set/config.yaml -j 64
```

#### Important
This workflow expects classified prokaryotic contigs to already exist, so run classify first.

---

### 4. Metagenome summarisation

**Workflow:** `summarise_metagenome_v01.smk`

This step currently:
- runs eggNOG annotation on prokaryotic contigs
- calculates multi-sample coverage summaries with CoverM

#### Dry run
```bash name=summarise_dry_run.sh
snakemake -s summarise_metagenome_v01.smk --use-conda --configfile test_set/config.yaml -n
```

#### Run
```bash name=summarise_run.sh
snakemake -s summarise_metagenome_v01.smk --use-conda --configfile test_set/config.yaml -j 64
```

#### Important
This workflow expects renamed prokaryotic contigs and assembly groups from earlier steps.

---

### 5. Bin summarisation (coming soon)

**Planned workflow:** `summarise_bins.smk`

This upcoming workflow will operate on recovered bins and is planned to include:
- dereplication using **Galah**
- taxonomy assignment with **GTDB-Tk**
- completeness / contamination estimation with **CheckM2**
- functional profiling with **METABOLIC-G**

Planned purpose:
- take binning outputs
- refine the MAG catalog
- annotate taxonomy, quality, and metabolic potential

This will become the recommended downstream step after `metagenome_binning_march2026.smk`.

---

## 🚀 Recommended full run order

```bash name=full_pipeline_tutorial.sh
snakemake -s metagenome_assemble_march2026.smk --use-conda --configfile test_set/config.yaml -j 64
snakemake -s metagenome_classify_march2026.smk --use-conda --configfile test_set/config.yaml -j 64
snakemake -s metagenome_binning_march2026.smk --use-conda --configfile test_set/config.yaml -j 64
snakemake -s summarise_metagenome_v01.smk --use-conda --configfile test_set/config.yaml -j 64
```

When `summarise_bins.smk` is added, it will be the next downstream step after binning.

---

## 📁 Example repository structure

```text name=repo_structure.txt
MACE_metagenomics/
├── README.md
├── metagenome_assemble_march2026.smk
├── metagenome_classify_march2026.smk
├── metagenome_binning_march2026.smk
├── summarise_metagenome_v01.smk
├── summarise_bins.smk                 # planned
├── envs/
│   ├── MACE_metagenomics.yml
│   ├── binette.yml
│   ├── binning_base.yml
│   ├── comebin.yml
│   ├── coverm.yml
│   ├── maxbin2.yml
│   ├── metabat2.yml
│   ├── semibin2.yml
│   ├── whokaryote.yml
│   └── ...
└── test_set/
    └── config.yaml
```

---

## 🧷 Practical notes

- Always use `--use-conda` unless you are managing environments yourself.
- Use absolute paths in config files when running on a server.
- A dry run with `-n` is strongly recommended before launching a full execution.
- If classify, binning, or summarise fail at parse time, check whether the upstream outputs already exist.
- Logs are written by individual rules and are the first place to inspect failures.

---

## 🛠️ Troubleshooting

- **Missing samples in metadata:** ensure the `Sample` column matches the Illumina file names.
- **Classification step fails early:** confirm that assembly has already produced `assembly_groups.tsv`.
- **Binning step cannot find prokaryotic contigs:** confirm that classification completed successfully.
- **Coverage-related errors:** check that the trimmed short reads exist in the configured `binning.short_reads_folder`.
- **Environment issues:** verify that the env YAML file referenced by each rule exists and can be solved by conda/mamba.

---

## 📣 Getting help

- Open a [GitHub issue](https://github.com/MACE-laboratory/MACE_metagenomics/issues) for bugs or feature requests.
- Check the workflow files directly for parameter expectations and output naming.

---

Happy assembling, classifying, binning, and summarising!

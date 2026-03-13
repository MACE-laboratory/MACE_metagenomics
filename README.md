# MACE Metagenomics: Snakemake Workflows

Welcome to the MACE Metagenomics repository!  
This project provides modular, scalable, and reproducible Snakemake workflows for metagenomic assembly, classification, and binning.

Workflows are designed for flexibility and ease-of-use, and each integrates conda environments for reproducibility.

---

## 📂 Available Workflows

| Workflow                                   | Purpose                                   |
|---------------------------------------------|-------------------------------------------|
| `metagenome_assemble_march2026.smk`        | Metagenome assembly from short/long reads |
| `metagenome_classify_march2026.smk`        | Eukaryote/prokaryote contig classification|
| `metagenome_binning_march2026.smk`         | Genome binning from metagenome assemblies |

All workflows use configs and environment YAMLs, making them modular and maintainable.

---

## 🔬 1. Metagenomic Assembly

**File:** `metagenome_assemble_march2026.smk`  
**Purpose:** Performs (optionally hybrid) assembly of metagenomic data from short (Illumina) and/or long (Nanopore) reads using MEGAHIT, MetaSPAdes, and others.  
**Inputs:**  
- Illumina FASTQ reads in a specified folder  
- (Optional) Nanopore FASTQ reads in a folder  
- Metadata table  
- Config YAML

**Output:**  
- Assembled contigs in output directory  
- Summary reports (QUAST, CheckM, etc.)

**Example config (`config_assemble.yml`):**
```yaml
outdir: "results"
illumina_folder: "data/short_reads"
nanopore_folder: "data/long_reads"   # Optional if only short reads
cpus: 8
trim_galore_threads: 4
metadata: "samples_metadata.tsv"
assembly_type: ["co", "sub", "single"]
long_reads_preprocessing:
  porechop: true
  chopper: false
```
| Parameter               | Explanation                                              |
|-------------------------|---------------------------------------------------------|
| `outdir`                | Output directory                                        |
| `illumina_folder`       | Folder (or glob) of Illumina reads                      |
| `nanopore_folder`       | Folder (or glob) of Nanopore reads (optional)           |
| `cpus`                  | Total CPUs to use                                       |
| `trim_galore_threads`   | Threads for qc/trimming                                 |
| `metadata`              | Tab or comma-separated sample metadata                  |
| `assembly_type`         | List: "co", "sub", "single" - which type(s) to run      |
| `long_reads_preprocessing` | Tools to apply to long reads (booleans for each)     |

**Run the workflow:**
```sh
snakemake --snakefile metagenome_assemble_march2026.smk --configfile config_assemble.yml --use-conda -j 8
```

---

## 🧬 2. Contig Classification (Eukaryote/Prokaryote)

**File:** `metagenome_classify_march2026.smk`  
**Purpose:** Classifies contigs in your assemblies as prokaryote or eukaryote using WhoKaryote.  
**Inputs:**  
- `outdir/assembly_groups.tsv` (from assembly workflow)
- All assemblies in `outdir/assemblies/`
- Config YAML

**Output:**  
- Table with classification results per contig

**Example config (`config_classify.yml`):**
```yaml
outdir: "results"
cpus: 4
```
| Parameter   | Explanation                           |
|-------------|--------------------------------------|
| `outdir`    | Directory with assemblies             |
| `cpus`      | Number of CPUs to use                 |

**Run the workflow:**
```sh
snakemake --snakefile metagenome_classify_march2026.smk --configfile config_classify.yml --use-conda -j 4
```

---

## 🧱 3. Genome Binning

**File:** `metagenome_binning_march2026.smk`  
**Purpose:** Performs binning on assemblies to recover MAGs using multiple binner tools (MetaBAT2, MaxBin2, Binette, SemiBin2, Rosella, etc.).  
**Inputs:**  
- Assembly FASTA file  
- Short reads directory  
- Config YAML

**Output:**  
- Binned genomes (MAGs) in output directory  
- Quality reports

**Example config (`config_binning.yml`):**
```yaml
outdir: "results"
threads: 8
assembly_path: "results/assemblies/my_sample.fa"
suffix: "mysample"
short_reads_folder: "results/trimmed_reads/short_reads"
```
| Parameter            | Explanation                                    |
|----------------------|------------------------------------------------|
| `outdir`             | Output directory for bins/results              |
| `threads`            | Number of threads to use                       |
| `assembly_path`      | FASTA file to bin                              |
| `suffix`             | Sample or condition name                       |
| `short_reads_folder` | Location of matched short reads                |

**Run the workflow:**
```sh
snakemake --snakefile metagenome_binning_march2026.smk --configfile config_binning.yml --use-conda -j 8
```

---

## 💡 Conda Environments

All required conda environment YAMLs are found in `workflow/envs/` and are used automatically.  
If you are running on a server/cluster, ensure conda is available and updated.

---

## 🧷 Example Folder Structure

```
MACE_metagenomics/
├── metagenome_assemble_march2026.smk
├── metagenome_binning_march2026.smk
├── metagenome_classify_march2026.smk
├── workflow/
│   └── envs/
│        ├── binette.yml
│        ├── maxbin2.yml
│        ├── metabat2.yml
│        ├── semibin2.yml
│        └── ... (other envs)
├── config_assemble.yml
├── config_classify.yml
├── config_binning.yml
└── README.md
```

---

## 🛠️ Tips & Troubleshooting

- Use `--use-conda` to auto-create environments per rule
- To run specific rules: `snakemake -s Snakefile -R <rulename>`
- For clusters, adapt the snakemake command line to your scheduler (see Snakemake docs for profiles)
- If something fails, consult the error log for missing files/typos in config
- Always use fresh YAMLs and fastq paths—absolute paths reduce confusion!

---

## 📣 Getting Help

- For problems/feature requests, open a [GitHub issue](https://github.com/MACE-laboratory/MACE_metagenomics/issues)
- For help understanding parameter choices, consult the comments at the top of each workflow file.

---

Happy assembling, classifying, and binning!

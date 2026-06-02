import os
import pandas as pd
from pathlib import Path
import glob

OUTDIR = config["outdir"]
SHORT_DIR = f"{OUTDIR}/trimmed_reads/short_reads"
LONG_DIR = f"{OUTDIR}/trimmed_reads/long_reads"
ASSEMBLY_DIR = f"{OUTDIR}/assemblies"

ILLUMINA = config["illumina_folder"]
NANOPORE = config.get("nanopore_folder")
IS_HYBRID = NANOPORE not in [None, "", "null"]

ASSEMBLY_TYPES = config["assembly_type"]

############################################
# Create directories
############################################

os.makedirs(SHORT_DIR, exist_ok=True)
os.makedirs(LONG_DIR, exist_ok=True)
os.makedirs(ASSEMBLY_DIR, exist_ok=True)

############################################
# Helper functions
############################################

def find_read_files(sample_name, reads_folder):
    """
    Find read files for a given sample name.
    Looks for files matching: sample_name_R1_001.f*q.gz and sample_name_R2_001.f*q.gz
    """
    r1_file = os.path.join(reads_folder, f"{sample_name}_R1_001.f*q.gz")
    r2_file = os.path.join(reads_folder, f"{sample_name}_R2_001.f*q.gz")
    
    r1_files = sorted(glob.glob(r1_file))
    r2_files = sorted(glob.glob(r2_file))
    
    return r1_files, r2_files

def list_short_reads():
    """List all R1 read files"""
    samples = glob.glob(os.path.join(ILLUMINA, "*_R1_001.f*q.gz"))
    samples = [s for s in samples if "control" not in s and "Nc" not in s]
    return samples

def sample_name_from_r1(path):
    """Extract sample name from R1 file path"""
    basename = os.path.basename(path)
    # Remove everything from _R1_001 onwards
    return basename.split("_R1_001")[0]

def validate_metadata_samples(metadata_file):
    """
    Validate that all samples in metadata file have corresponding raw read files
    """
    df = pd.read_csv(metadata_file, sep=None, engine="python")
    
    # Get available raw samples
    raw_samples_list = list_short_reads()
    raw_samples = set(sample_name_from_r1(s) for s in raw_samples_list)
    
    # Get samples from metadata
    metadata_samples = set(df["Sample"].tolist())
    
    # Check for missing samples
    missing = metadata_samples - raw_samples
    
    if missing:
        raise ValueError(
            f"The following samples in metadata file '{metadata_file}' do not have "
            f"corresponding raw read files in '{ILLUMINA}':\n{missing}\n\n"
            f"Available samples: {raw_samples}"
        )
    
    # Validate that each sample has both R1 and R2 reads
    for sample in metadata_samples:
        r1_files, r2_files = find_read_files(sample, ILLUMINA)
        
        if not r1_files:
            raise ValueError(
                f"Sample '{sample}' from metadata has no R1 read files matching "
                f"pattern '{sample}_R1_001.f*q.gz' in '{ILLUMINA}'"
            )
        
        if not r2_files:
            raise ValueError(
                f"Sample '{sample}' from metadata has no R2 read files matching "
                f"pattern '{sample}_R2_001.f*q.gz' in '{ILLUMINA}'"
            )
        
        if len(r1_files) != len(r2_files):
            raise ValueError(
                f"Sample '{sample}' has mismatched number of R1 and R2 files:\n"
                f"  R1 files ({len(r1_files)}): {r1_files}\n"
                f"  R2 files ({len(r2_files)}): {r2_files}"
            )
    
    return df

SHORT_READS = list_short_reads()
SAMPLES = sorted(list(set(sample_name_from_r1(r) for r in SHORT_READS)))

# Validate metadata on workflow start
METADATA_DF = validate_metadata_samples(config["metadata"])

############################################
# Rule: all
############################################

def get_final_outputs(wildcards):
    ck = checkpoints.parse_metadata.get()
    groups = ck.output.groups
    df = pd.read_csv(groups, sep="\t")
    outputs = df["assembly_path"].tolist()
    # Add the tools versions file
    outputs.append(f"{ASSEMBLY_DIR}/tools_versions.yaml")
    return outputs

rule all:
    input:
        get_final_outputs

############################################
# Trim short reads
############################################

rule trim_galore:
    input:
        r1_files=lambda wildcards: find_read_files(wildcards.sample, ILLUMINA)[0],
        r2_files=lambda wildcards: find_read_files(wildcards.sample, ILLUMINA)[1]
    output:
        r1=f"{SHORT_DIR}/{{sample}}_R1_001_val_1.fq.gz",
        r2=f"{SHORT_DIR}/{{sample}}_R2_001_val_2.fq.gz"
    threads:
        config["trim_galore_threads"]
    log:
        "logs/trim_galore_{sample}.log"
    conda:
        "envs/MACE_metagenomics.yml"
    shell:
        """
        trim_galore --fastqc --paired -j {threads} \
        -o {SHORT_DIR} {input.r1_files} {input.r2_files} 2>&1 | tee {log}
        """

############################################
# Mark short reads trimming complete
############################################

rule trim_short_reads_done:
    input:
        expand(f"{SHORT_DIR}/{{sample}}_R1_001_val_1.fq.gz", sample=SAMPLES),
        expand(f"{SHORT_DIR}/{{sample}}_R2_001_val_2.fq.gz", sample=SAMPLES)
    output:
        f"{OUTDIR}/trimmed_reads/short_reads.done"
    conda:
        "envs/MACE_metagenomics.yml"
    run:
        # Verify that all trimmed files exist and are valid
        missing_files = []
        for sample in SAMPLES:
            r1 = f"{SHORT_DIR}/{sample}_R1_001_val_1.fq.gz"
            r2 = f"{SHORT_DIR}/{sample}_R2_001_val_2.fq.gz"
            
            if not os.path.exists(r1) or not os.path.exists(r2):
                missing_files.append((r1, r2))
        
        if missing_files:
            raise ValueError(
                f"The following trimmed read files are missing:\n" +
                "\n".join([f"  {r1}\n  {r2}" for r1, r2 in missing_files])
            )
        
        shell("touch {output}")

############################################
# Long reads preprocessing (flexible)
############################################

rule preprocess_long_reads:
    input:
        reads=lambda wc: wc.read
    output:
        f"{LONG_DIR}/{{sample}}_chopped.fastq.gz"
    threads:
        config["cpus"]
    log:
        "logs/preprocess_long_reads_{sample}.log"
    conda:
        "envs/MACE_metagenomics.yml"
    run:
        inp = input.reads
        out = output[0]

        if config["long_reads_preprocessing"]["porechop"]:
            pore = out.replace("_chopped.fastq.gz", "_porechopped.fastq.gz")
            shell(f"porechop --threads {threads} -i {inp} -o {pore} 2>&1 | tee {log}")
            inp = pore

        if config["long_reads_preprocessing"]["chopper"]:
            shell(f"gunzip -c {inp} | chopper -q 12 --threads {threads} | gzip > {out} 2>&1 | tee {log}")
        else:
            shell(f"cp {inp} {out}")

############################################
# Mark long reads preprocessing complete
############################################

rule trim_long_reads_done:
    input:
        expand(f"{LONG_DIR}/{{sample}}_chopped.fastq.gz", sample=SAMPLES) if IS_HYBRID else []
    output:
        f"{OUTDIR}/trimmed_reads/long_reads.done"
    conda:
        "envs/MACE_metagenomics.yml"
    run:
        if IS_HYBRID:
            # Verify that all long read files exist and are valid
            missing_files = []
            for sample in SAMPLES:
                long_read = f"{LONG_DIR}/{sample}_chopped.fastq.gz"
                
                if not os.path.exists(long_read):
                    missing_files.append(long_read)
            
            if missing_files:
                raise ValueError(
                    f"The following long read files are missing:\n" +
                    "\n".join([f"  {f}" for f in missing_files])
                )
        
        shell("touch {output}")

############################################
# Checkpoint: parse metadata and decide assemblies
############################################

checkpoint parse_metadata:
    input:
        metadata=config["metadata"]
    output:
        groups=f"{OUTDIR}/assembly_groups.tsv"
    conda:
        "envs/MACE_metagenomics.yml"
    run:
        df = pd.read_csv(input.metadata, sep=None, engine="python")

        results = []

        if "co" in ASSEMBLY_TYPES:
            results.append({
                "assembly": "coassembly",
                "samples": ",".join(df["Sample"].tolist())
            })

        if "sub" in ASSEMBLY_TYPES:
            for g, subdf in df.groupby("Group"):
                results.append({
                    "assembly": f"sub_assembly_{g}",
                    "samples": ",".join(subdf["Sample"].tolist())
                })

        if "single" in ASSEMBLY_TYPES:
            for s in df["Sample"]:
                results.append({
                    "assembly": f"single_assembly_{s}",
                    "samples": s
                })

        out = pd.DataFrame(results)

        paths = []
        for a in out["assembly"]:
            if IS_HYBRID:
                paths.append(f"{ASSEMBLY_DIR}/{a}/contigs.fasta")
            else:
                paths.append(f"{ASSEMBLY_DIR}/{a}/final.contigs.fa")

        out["assembly_path"] = paths
        out.to_csv(output.groups, sep="\t", index=False)

############################################
# Helper function to get samples for an assembly
############################################

def get_samples_for_assembly(assembly_name):
    """Extract sample names from assembly_groups.tsv for a given assembly"""
    try:
        ck = checkpoints.parse_metadata.get()
        groups = ck.output.groups
        df = pd.read_csv(groups, sep="\t")
        row = df[df["assembly"] == assembly_name]
        if len(row) > 0:
            return row.iloc[0]["samples"].split(",")
        return []
    except:
        return []

############################################
# Short-read assembly (MEGAHIT)
############################################

rule megahit:
    input:
        short_reads_done=f"{OUTDIR}/trimmed_reads/short_reads.done",
        groups=lambda wildcards: checkpoints.parse_metadata.get().output.groups
    output:
        contigs=f"{ASSEMBLY_DIR}/{{assembly}}/final.contigs.fa"
    threads:
        config["cpus"]
    log:
        "logs/megahit_{assembly}.log"
    conda:
        "envs/MACE_metagenomics.yml"
    run:
        samples = get_samples_for_assembly(wildcards.assembly)
        
        # Validate that trimmed reads exist for all samples in this assembly
        missing_reads = []
        for sample in samples:
            r1 = f"{SHORT_DIR}/{sample}_R1_001_val_1.fq.gz"
            r2 = f"{SHORT_DIR}/{sample}_R2_001_val_2.fq.gz"
            
            if not os.path.exists(r1):
                missing_reads.append(r1)
            if not os.path.exists(r2):
                missing_reads.append(r2)
        
        if missing_reads:
            raise ValueError(
                f"Cannot run megahit for assembly '{wildcards.assembly}'. "
                f"The following trimmed read files are missing:\n" +
                "\n".join([f"  {f}" for f in missing_reads])
            )
        
        fwd_files = ",".join([f"{SHORT_DIR}/{s}_R1_001_val_1.fq.gz" for s in samples])
        rev_files = ",".join([f"{SHORT_DIR}/{s}_R2_001_val_2.fq.gz" for s in samples])
        out_dir = f"{ASSEMBLY_DIR}/{wildcards.assembly}"

        shell(f"rm -rf {out_dir} && megahit -m 0.95 --presets meta-large --k-min 27 --k-max 127 --k-step 10 --min-contig-len 1000 -t {threads} -1 {fwd_files} -2 {rev_files} -o {out_dir} 2>&1 | tee {log}")

############################################
# Hybrid assembly (MetaSPAdes)
############################################

rule metaspades_hybrid:
    input:
        short_reads_done=f"{OUTDIR}/trimmed_reads/short_reads.done",
        long_reads_done=f"{OUTDIR}/trimmed_reads/long_reads.done",
        groups=lambda wildcards: checkpoints.parse_metadata.get().output.groups
    output:
        contigs=f"{ASSEMBLY_DIR}/{{assembly}}/contigs.fasta"
    threads:
        config["cpus"]
    log:
        "logs/metaspades_{assembly}.log"
    conda:
        "envs/MACE_metagenomics.yml"
    run:
        samples = get_samples_for_assembly(wildcards.assembly)
        
        # Validate that trimmed reads exist for all samples in this assembly
        missing_reads = []
        for sample in samples:
            r1 = f"{SHORT_DIR}/{sample}_R1_001_val_1.fq.gz"
            r2 = f"{SHORT_DIR}/{sample}_R2_001_val_2.fq.gz"
            long_r = f"{LONG_DIR}/{sample}_chopped.fastq.gz"
            
            if not os.path.exists(r1):
                missing_reads.append(r1)
            if not os.path.exists(r2):
                missing_reads.append(r2)
            if not os.path.exists(long_r):
                missing_reads.append(long_r)
        
        if missing_reads:
            raise ValueError(
                f"Cannot run metaspades for assembly '{wildcards.assembly}'. "
                f"The following trimmed read files are missing:\n" +
                "\n".join([f"  {f}" for f in missing_reads])
            )
        
        assembly_dir = f"{ASSEMBLY_DIR}/{wildcards.assembly}"
        fwd_files = " ".join([f"{SHORT_DIR}/{s}_R1_001_val_1.fq.gz" for s in samples])
        rev_files = " ".join([f"{SHORT_DIR}/{s}_R2_001_val_2.fq.gz" for s in samples])
        long_files = " ".join([f"{LONG_DIR}/{s}_chopped.fastq.gz" for s in samples])

        shell(f"rm -rf {assembly_dir} && cat {fwd_files} > /tmp/all_forward_{wildcards.assembly}.fastq.gz && cat {rev_files} > /tmp/all_reverse_{wildcards.assembly}.fastq.gz && cat {long_files} > /tmp/all_long_{wildcards.assembly}.fastq.gz && spades.py -t {threads} -m 950 -1 /tmp/all_forward_{wildcards.assembly}.fastq.gz -2 /tmp/all_reverse_{wildcards.assembly}.fastq.gz --nanopore /tmp/all_long_{wildcards.assembly}.fastq.gz -o {assembly_dir} 2>&1 | tee {log}")

############################################
# Generate tools versions file
############################################

rule generate_tools_versions:
    input:
        lambda wildcards: checkpoints.parse_metadata.get().output.groups
    output:
        versions=f"{ASSEMBLY_DIR}/tools_versions.yaml"
    conda:
        "envs/MACE_metagenomics.yml"
    shell:
        """
        cat > {output.versions} << 'EOFHEADER'
tools_versions:
  trim_galore:
    command: "trim_galore --version"
    description: "Quality control and adapter trimming for reads"
  
  megahit:
    command: "megahit --version"
    description: "Short-read assembly tool"
  
  spades:
    command: "spades.py --version"
    description: "Hybrid assembly tool (MetaSPAdes)"
  
  porechop:
    command: "porechop --version"
    description: "Long-read adapter trimming tool"
    optional: true
  
  chopper:
    command: "chopper --version"
    description: "Long-read quality filtering tool"
    optional: true

version_generation_date: "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
execution_command: "snakemake metagenome_assembly_march2026.smk --use-conda"

EOFHEADER

        echo "" >> {output.versions}
        echo "# Actual versions captured at execution time:" >> {output.versions}
        echo "actual_versions:" >> {output.versions}
        
        trim_galore_version=$(trim_galore --version 2>&1 | head -1) && echo "  trim_galore: $trim_galore_version" >> {output.versions} || echo "  trim_galore: version unknown" >> {output.versions}
        megahit_version=$(megahit --version 2>&1 | head -1) && echo "  megahit: $megahit_version" >> {output.versions} || echo "  megahit: version unknown" >> {output.versions}
        spades_version=$(spades.py --version 2>&1 | head -1) && echo "  spades: $spades_version" >> {output.versions} || echo "  spades: version unknown" >> {output.versions}
        
        if command -v porechop &> /dev/null; then
            porechop_version=$(porechop --version 2>&1 | head -1) && echo "  porechop: $porechop_version" >> {output.versions} || echo "  porechop: version unknown" >> {output.versions}
        fi
        
        if command -v chopper &> /dev/null; then
            chopper_version=$(chopper --version 2>&1 | head -1) && echo "  chopper: $chopper_version" >> {output.versions} || echo "  chopper: version unknown" >> {output.versions}
        fi
        """
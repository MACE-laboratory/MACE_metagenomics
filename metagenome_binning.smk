import os
import glob
import pandas as pd

OUTDIR = config["outdir"]
RESOURCES = config["resources"]
BINNING_CONFIG = config["binning"]

threads_n = RESOURCES["threads"]
short_reads_folder = f"{OUTDIR}/trimmed_reads/short_reads"

if "suffixes" in BINNING_CONFIG and BINNING_CONFIG["suffixes"]:
    suffixes = BINNING_CONFIG["suffixes"]
else:
    groups = pd.read_csv(f"{OUTDIR}/assembly_groups.tsv", sep="\t", dtype=str)
    suffixes = groups["assembly"].tolist()

fastq_files = glob.glob(f"{short_reads_folder}/*R1_001_val_1.fq.gz")
samples = [os.path.basename(f).replace("_R1_001_val_1.fq.gz", "") for f in fastq_files]

def get_assembly(wildcards):
    return f"{OUTDIR}/assemblies/{wildcards.suffix}/prokaryotic_contigs.fasta"

def get_bam_file_path(wildcards):
    checkpoint_output = checkpoints.run_coverm.get(suffix=wildcards.suffix).output.bam_files
    assembly_name = "prokaryotic_contigs.renamed.filtered.fa"
    bam_filename = f"{assembly_name}.{wildcards.sample}_R1_001_val_1.fq.gz.bam"
    return os.path.join(checkpoint_output, bam_filename)

rule all:
    input:
        expand(
            f"{OUTDIR}/binning/{{suffix}}/binette/final_bins_quality_reports.tsv",
            suffix=suffixes
        )

rule create_binning_folder:
    output:
        done = touch(f"{OUTDIR}/binning/{{suffix}}/.binning_folder.done")
    shell:
        """
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}/coverm
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}/maxbin2
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}/semibin2
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}/metabat2
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}/binette
        mkdir -p {OUTDIR}/binning/{{wildcards.suffix}}/comebin
        touch {output.done}
        """

rule filter_contigs:
    input:
        get_assembly
    output:
        f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.filtered.fa"
    conda:
        "envs/binning_base.yml"
    shell:
        "seqkit seq -j 10 --remove-gaps -m 1500 {input} > {output}"

rule rename_contigs:
    input:
        f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.filtered.fa"
    output:
        f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"
    conda:
        "envs/binning_base.yml"
    shell:
        "perl -pe 's/ .*//g' {input} > {output}"

checkpoint run_coverm:
    input:
        assembly = f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        binning_folders = f"{OUTDIR}/binning/{{suffix}}/.binning_folder.done"
    output:
        coverm_metabat = f"{OUTDIR}/binning/{{suffix}}/coverm/coverm_metabat.tsv",
        bam_files = directory(f"{OUTDIR}/binning/{{suffix}}/coverm/bam_files/bam")
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        """
        coverm contig \
        -p bwa-mem2 \
        -r {input.assembly} \
        -m metabat \
        -t {threads} \
        -1 {short_reads_folder}/*val_1.fq.gz \
        -2 {short_reads_folder}/*val_2.fq.gz \
        -o {output.coverm_metabat} \
        --bam-file-cache-directory {output.bam_files}
        """

rule coverage_maxbin:
    input:
        bam=get_bam_file_path,
        binning_folders = f"{OUTDIR}/binning/{{suffix}}/.binning_folder.done"
    output:
        f"{OUTDIR}/binning/{{suffix}}/coverm/bam_files/trimmed_means/{{sample}}.trimmed_mean"
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        "coverm contig -t {threads} -m trimmed_mean -b {input.bam} > {output}"

rule list_trimmed_files:
    input:
        lambda wildcards: expand(
            f"{OUTDIR}/binning/{{suffix}}/coverm/bam_files/trimmed_means/{{sample}}.trimmed_mean",
            sample=samples,
            suffix=wildcards.suffix
        )
    output:
        f"{OUTDIR}/binning/{{suffix}}/coverm/trimmed_files_list.txt"
    shell:
        "printf '%s\\n' {input} > {output}"

rule run_maxbin2:
    input:
        assembly = f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        coverm = f"{OUTDIR}/binning/{{suffix}}/coverm/trimmed_files_list.txt"
    output:
        f"{OUTDIR}/binning/{{suffix}}/maxbin2/maxbin2_run.summary"
    params:
        out_prefix = f"{OUTDIR}/binning/{{suffix}}/maxbin2/maxbin2_run"
    conda:
        "envs/maxbin2.yml"
    threads: threads_n
    shell:
        "run_MaxBin.pl -contig {input.assembly} -abund_list {input.coverm} -out {params.out_prefix} -min_contig_length 1500 -thread {threads}"

rule gzip_assembly:
    input:
        f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"
    output:
        f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa.gz"
    shell:
        "gzip -c {input} > {output}"

rule run_metabat2:
    input:
        assembly = f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa.gz",
        coverm = f"{OUTDIR}/binning/{{suffix}}/coverm/coverm_metabat.tsv"
    params:
        out_path = f"{OUTDIR}/binning/{{suffix}}/metabat2/bins_metabat"
    output:
        f"{OUTDIR}/binning/{{suffix}}/metabat2/bins_metabat.1.fa"
    conda:
        "envs/metabat2.yml"
    threads: threads_n
    shell:
        "metabat2 -i {input.assembly} -o {params.out_path} -a {input.coverm} -t {threads} --minContig 1500 --seed 23"

rule rename_metabat2:
    input:
        f"{OUTDIR}/binning/{{suffix}}/metabat2/bins_metabat.1.fa"
    params:
        dir = f"{OUTDIR}/binning/{{suffix}}/metabat2/"
    output:
        f"{OUTDIR}/binning/{{suffix}}/metabat2/bins_metabat_1.fa"
    shell:
        "cp {input} {output} && rename 's/bins_metabat\\./bins_metabat_/g' {params.dir}*.fa"

rule run_comebin:
    input:
        assembly = f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        bam_dir = f"{OUTDIR}/binning/{{suffix}}/coverm/bam_files/bam"
    output:
        f"{OUTDIR}/binning/{{suffix}}/comebin/comebin_res/comebin_res.tsv"
    params:
        outdir = f"{OUTDIR}/binning/{{suffix}}/comebin"
    conda:
        "envs/comebin.yml"
    threads: threads_n
    shell:
        "CUDA_VISIBLE_DEVICES='' && run_comebin.sh -a {input.assembly} -o {params.outdir} -p {input.bam_dir} -t {threads}"

rule run_binette:
    input:
        assembly = f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        maxbin2 = f"{OUTDIR}/binning/{{suffix}}/maxbin2/maxbin2_run.summary",
        metabat2 = f"{OUTDIR}/binning/{{suffix}}/metabat2/bins_metabat_1.fa",
        comebin = f"{OUTDIR}/binning/{{suffix}}/comebin/comebin_res/comebin_res.tsv"
    params:
        metabat2 = f"{OUTDIR}/binning/{{suffix}}/metabat2",
        maxbin2 = f"{OUTDIR}/binning/{{suffix}}/maxbin2",
        comebin = f"{OUTDIR}/binning/{{suffix}}/comebin/comebin_res/comebin_res_bins",
        out_prefix = f"{OUTDIR}/binning/{{suffix}}/binette"
    output:
        f"{OUTDIR}/binning/{{suffix}}/binette/final_bins_quality_reports.tsv"
    threads: threads_n
    conda:
        "envs/binette.yml"
    shell:
        """
        export CHECKM2DB='/data/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd'
        binette --verbose \
        --bin_dirs {params.metabat2} {params.maxbin2} {params.comebin} \
        -c {input.assembly} \
        -t {threads} \
        -o {params.out_prefix}
        """

import os
import glob

# Configuration
assembly_path = config["assembly_path"]
assembly, fasta_ext = os.path.splitext(assembly_path)
assembly_renamed_filtered = assembly + '.renamed.filtered.fa'
outdir = config["outdir"]
suffix = config["suffix"]
threads_n = config["threads"]
short_reads_folder = config["short_reads_folder"]

fastq_files = glob.glob(f"{short_reads_folder}/*R1_001_val_1.fq.gz")
samples = [f.split("/")[-1].replace("_R1_001_val_1.fq.gz", "") for f in fastq_files]
print(samples)


def get_bam_file_paths(wildcards):
    checkpoint_output = checkpoints.run_coverm.get(**wildcards).output[1]
    # Construct specific BAM path based on wildcards.sample
    return os.path.join(checkpoint_output, f"{os.path.basename(assembly_path)}.{wildcards.sample}_R1_001_val_1.fq.gz.bam")


rule all:
    input:
        #f"{outdir}/binning/{suffix}/maxbin2/maxbin2_run.summary",
        #f"{outdir}/binning/{suffix}/metabat2/bins_metabat/bins_metabat.1.fa",
        #f"{outdir}/binning/{suffix}/rosella/rosella_bin_1.fna",
        f"{outdir}/binning/{suffix}/binette/final_bins_quality_reports.tsv"

rule create_binning_folder: 
    output:
        directory(f"{outdir}/binning/{suffix}")
    conda:
        f"envs/binning_base.yml"
    shell:
        "mkdir -p {output} {output}/coverm {output}/maxbin2 {output}/semibin2 {output}/metabat2 {output}/rosella {output}/binette"

rule filter_contigs:
    input:
        f"{assembly_path}"
    output:
        f"{assembly}.filtered.fa"
    conda:
        "envs/binning_base.yml"
    shell:
        "seqkit seq -j 10 --remove-gaps -o {output} -m 1500 {input}"

rule rename_contigs:
    input:
        f"{assembly}.filtered.fa"    
    output:
        f"{assembly}.renamed.filtered.fa"
    conda:
        "envs/binning_base.yml"
    shell:
        "perl -pe 's/ .*//g' {input} > {output}"

checkpoint run_coverm:
    input:
        assembly=f'{assembly}.renamed.filtered.fa'
    output:
        coverm_metabat=f"{outdir}/binning/{suffix}/coverm/coverm_metabat.tsv", 
        bam_files=directory(f"{outdir}/binning/{suffix}/coverm/bam_files/bam")
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        f"coverm contig -p bwa-mem2 -r {input.assembly} -m metabat -t {threads_n} -1 {short_reads_folder}/*val_1.fq.gz -2 {short_reads_folder}/*val_2.fq.gz -o {output.coverm_metabat} --bam-file-cache-directory {output.bam_files} --discard-unmapped"

rule coverage_maxbin:
    input:
        bam = get_bam_file_paths
    output:
        f"{outdir}/binning/{suffix}/coverm/bam_files/trimmed_means/{{sample}}.trimmed_mean"
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        "coverm contig -t {threads} -m trimmed_mean -b {input.bam} > {output}" 

rule list_trimmed_files:
    input:
        expand(f"{outdir}/binning/{suffix}/coverm/bam_files/trimmed_means/{{sample}}.trimmed_mean", sample=samples)
    output:
        f"{outdir}/binning/{suffix}/coverm/trimmed_files_list.txt"
    shell:
        "printf '%s\n' {input} > {output}"

rule run_maxbin2:
    input:
        assembly=f'{assembly}.renamed.filtered.fa',
        coverm=f"{outdir}/binning/{suffix}/coverm/trimmed_files_list.txt"
    output:
        f"{outdir}/binning/{{suffix}}/maxbin2/maxbin2_run.summary"
    params:
        out_prefix=lambda wildcards: f"{outdir}/binning/{wildcards.suffix}/maxbin2/maxbin2_run"
    conda:
        "envs/maxbin2.yml"
    threads: threads_n
    shell:
        "run_MaxBin.pl -contig {input.assembly} -abund_list {input.coverm} -out {params.out_prefix} -min_contig_length 1500 -thread {threads}"

rule run_rosella:
    input:
        assembly=f'{assembly}.renamed.filtered.fa',
        coverm=f'{outdir}/binning/{suffix}/coverm/coverm_metabat.tsv'
    params:
        outdir=directory(f"{outdir}/binning/{suffix}/rosella")
    output:
        f"{outdir}/binning/{suffix}/rosella/rosella_bin_unbinned.fna"
    conda:
        "envs/rosella.yml"
    threads: threads_n
    shell:
        f"rosella recover -C {input.coverm} -r {input.assembly} --output-directory {params.outdir} --threads {threads}"

rule create_rosella_dir:
    input:
        f"{outdir}/binning/{suffix}/rosella/rosella_bin_unbinned.fna"
    params:
        out_dir=f"{outdir}/binning/{suffix}/rosella/rosella_bins",
        to_copy=f"{outdir}/binning/{suffix}/rosella/*.fna"
    output:
        directory(f"{outdir}/binning/{suffix}/rosella/rosella_bins/"),
        touch(f"{outdir}/binning/{suffix}/rosella/rosella_bins/.done")
    shell:
        "cp {params.to_copy} {params.out_dir}/. && rm {params.out_dir}/*unbinned*.fna"

rule gzip_assembly:
    input:
        assembly=f'{assembly}.renamed.filtered.fa'
    output:
        assembly=f'{assembly}.renamed.filtered.fa.gz'
    shell:
        "gzip -c {input.assembly} > {output}"

rule run_metabat2:
    input:
        assembly=f'{assembly}.renamed.filtered.fa.gz',
        rosella_done=f"{outdir}/binning/{suffix}/rosella/rosella_bins/.done"
    params:
        out_path=f"{outdir}/binning/{suffix}/metabat2/bins_metabat"
    output:
        f"{outdir}/binning/{suffix}/metabat2/bins_metabat.1.fa"
    conda:
        "envs/metabat2.yml"
    threads: threads_n
    shell:
        f"metabat2 -i {input.assembly} -o {params.out_path} -a {outdir}/binning/{{suffix}}/coverm/coverm_metabat.tsv -t {threads_n} --minContig 1500 --seed 23"

rule rename_metabat2:
    input:
        f"{outdir}/binning/{suffix}/metabat2/bins_metabat.1.fa"
    params:
        dir=f"{outdir}/binning/{suffix}/metabat2/"
    output:
        f"{outdir}/binning/{suffix}/metabat2/bins_metabat_1.fa"
    shell:
        "rename 's/bins_metabat\\./bins_metabat_/g' {params.dir}*.fa"

rule run_binette:
    input:
        rosella=f"{outdir}/binning/{suffix}/rosella/rosella_bins/.done",
        maxbin2=f"{outdir}/binning/{suffix}/maxbin2/maxbin2_run.summary",
        metabat2=f"{outdir}/binning/{suffix}/metabat2/bins_metabat_1.fa",
        assembly=f"{assembly}.renamed.filtered.fa"
    params:
        rosella=f"{outdir}/binning/{suffix}/rosella/rosella_bins",
        metabat2=f"{outdir}/binning/{suffix}/metabat2",
        maxbin2=f"{outdir}/binning/{suffix}/maxbin2",
        out_prefix=f"{outdir}/binning/{suffix}/binette"
    output:
        f"{outdir}/binning/{suffix}/binette/final_bins_quality_reports.tsv"
    threads: threads_n
    conda:
        "envs/binette.yml"
    shell:
        f"export CHECKM2DB='/data/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd' && binette --verbose --bin_dirs {params.rosella} {params.metabat2} {params.maxbin2} -c {input.assembly} -t {threads_n} -o {params.out_prefix}"

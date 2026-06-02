import os
import glob

outdir = config["outdir"]
suffixes = config["suffixes"]
threads_n = config["threads"]
short_reads_folder = config["short_reads_folder"]

fastq_files = glob.glob(f"{short_reads_folder}/*R1_001_val_1.fq.gz")
samples = [os.path.basename(f).replace("_R1_001_val_1.fq.gz", "") for f in fastq_files]

def get_assembly(wildcards):
    return f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.fasta"

def get_filtered_assembly(wildcards):
    return f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.filtered.fa"

def get_renamed_assembly(wildcards):
    return f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"

def get_gzipped_assembly(wildcards):
    return f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa.gz"

def get_bam_file_path(wildcards):
    checkpoint_output = checkpoints.run_coverm.get(suffix=wildcards.suffix).output.bam_files
    assembly_name = "prokaryotic_contigs.renamed.filtered.fa"
    bam_filename = f"{assembly_name}.{wildcards.sample}_R1_001_val_1.fq.gz.bam"
    return os.path.join(checkpoint_output, bam_filename)

rule all:
    input:
        expand(
            f"{outdir}/binning/{{suffix}}/binette/final_bins_quality_reports.tsv",
            suffix=suffixes
        )

rule create_binning_folder:
    output:
        done = touch(f"{outdir}/binning/{{suffix}}/.binning_folder.done")
    shell:
        """
        mkdir -p {outdir}/binning/{{suffix}}
        mkdir -p {outdir}/binning/{{suffix}}/coverm
        mkdir -p {outdir}/binning/{{suffix}}/maxbin2
        mkdir -p {outdir}/binning/{{suffix}}/semibin2
        mkdir -p {outdir}/binning/{{suffix}}/metabat2
#        mkdir -p {outdir}/binning/{{suffix}}/rosella
        mkdir -p {outdir}/binning/{{suffix}}/binette
#        mkdir -p {outdir}/binning/{{suffix}}/comebin_res
        touch {output.done}
        """

rule filter_contigs:
    input:
        get_assembly
    output:
        f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.filtered.fa"
    conda:
        "envs/binning_base.yml"
    shell:
        "seqkit seq -j 10 --remove-gaps -m 1500 {input} > {output}"

rule rename_contigs:
    input:
        f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.filtered.fa"
    output:
        f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"
    conda:
        "envs/binning_base.yml"
    shell:
        "perl -pe 's/ .*//g' {input} > {output}"

checkpoint run_coverm:
    input:
        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        binning_folders = f"{outdir}/binning/{{suffix}}/.binning_folder.done"
    output:
        coverm_metabat = f"{outdir}/binning/{{suffix}}/coverm/coverm_metabat.tsv",
        bam_files = directory(f"{outdir}/binning/{{suffix}}/coverm/bam_files/bam")
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
        binning_folders = f"{outdir}/binning/{{suffix}}/.binning_folder.done"
    output:
        f"{outdir}/binning/{{suffix}}/coverm/bam_files/trimmed_means/{{sample}}.trimmed_mean"
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        "coverm contig -t {threads} -m trimmed_mean -b {input.bam} > {output}"

rule list_trimmed_files:
    input:
        lambda wildcards: expand(
            f"{outdir}/binning/{{suffix}}/coverm/bam_files/trimmed_means/{{sample}}.trimmed_mean",
            sample=samples,
            suffix=wildcards.suffix
        )
    output:
        f"{outdir}/binning/{{suffix}}/coverm/trimmed_files_list.txt"
    shell:
        "printf '%s\\n' {input} > {output}"

rule run_maxbin2:
    input:
        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        coverm = f"{outdir}/binning/{{suffix}}/coverm/trimmed_files_list.txt"
    output:
        f"{outdir}/binning/{{suffix}}/maxbin2/maxbin2_run.summary"
    params:
        out_prefix = f"{outdir}/binning/{{suffix}}/maxbin2/maxbin2_run"
    conda:
        "envs/maxbin2.yml"
    threads: threads_n
    shell:
        "run_MaxBin.pl -contig {input.assembly} -abund_list {input.coverm} -out {params.out_prefix} -min_contig_length 1500 -thread {threads}"

#rule run_rosella:
#    input:
#        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
#        coverm = f"{outdir}/binning/{{suffix}}/coverm/coverm_metabat.tsv"
#    params:
#        outdir = f"{outdir}/binning/{{suffix}}/rosella"
#    output:
#        f"{outdir}/binning/{{suffix}}/rosella/rosella_bin_unbinned.fna"
#    conda:
#        "test333"
#        #"envs/rosella.yml"
#    threads: threads_n
#    shell:
#        "rosella recover -C {input.coverm} -r {input.assembly} --output-directory {params.outdir} --threads {threads}"

#rule create_rosella_dir:
#    input:
#        f"{outdir}/binning/{{suffix}}/rosella/rosella_bin_unbinned.fna"
#    params:
#        out_dir = f"{outdir}/binning/{{suffix}}/rosella/rosella_bins",
#        to_copy = f"{outdir}/binning/{{suffix}}/rosella/*.fna"
#    output:
#        bins = directory(f"{outdir}/binning/{{suffix}}/rosella/rosella_bins/"),
#        done = touch(f"{outdir}/binning/{{suffix}}/rosella/rosella_bins/.done")
#    shell:
#        "mkdir -p {params.out_dir} && cp {params.to_copy} {params.out_dir}/. && rm {params.out_dir}/*unbinned*.fna && touch {output.done}"

rule gzip_assembly:
    input:
        f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"
    output:
        f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa.gz"
    shell:
        "gzip -c {input} > {output}"

rule run_metabat2:
    input:
        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa.gz",
        coverm = f"{outdir}/binning/{{suffix}}/coverm/coverm_metabat.tsv"#,
 #       rosella_done = f"{outdir}/binning/{{suffix}}/rosella/rosella_bins/.done"
    params:
        out_path = f"{outdir}/binning/{{suffix}}/metabat2/bins_metabat"
    output:
        f"{outdir}/binning/{{suffix}}/metabat2/bins_metabat.1.fa"
    conda:
        "envs/metabat2.yml"
    threads: threads_n
    shell:
        "metabat2 -i {input.assembly} -o {params.out_path} -a {input.coverm} -t {threads} --minContig 1500 --seed 23"

rule rename_metabat2:
    input:
        f"{outdir}/binning/{{suffix}}/metabat2/bins_metabat.1.fa"
    params:
        dir = f"{outdir}/binning/{{suffix}}/metabat2/"
    output:
        f"{outdir}/binning/{{suffix}}/metabat2/bins_metabat_1.fa"
    shell:
        "cp {input} {output} && rename 's/bins_metabat\\./bins_metabat_/g' {params.dir}*.fa"

rule run_comebin:
    input:
        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        bam_dir = f"{outdir}/binning/{{suffix}}/coverm/bam_files/bam"
    output:
        f"{outdir}/binning/{{suffix}}/comebin/comebin_res/comebin_res.tsv"
    params:
        outdir = f"{outdir}/binning/{{suffix}}/comebin"
    conda:
        "comebin_env"
    threads: threads_n
    shell:
        "CUDA_VISIBLE_DEVICES='' && run_comebin.sh -a {input.assembly} -o {params.outdir} -p {input.bam_dir} -t {threads}"

rule run_binette:
    input:
        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa",
        maxbin2 = f"{outdir}/binning/{{suffix}}/maxbin2/maxbin2_run.summary",
        metabat2 = f"{outdir}/binning/{{suffix}}/metabat2/bins_metabat_1.fa",
        comebin = f"{outdir}/binning/{{suffix}}/comebin/comebin_res/comebin_res.tsv"
    params:
        metabat2 = f"{outdir}/binning/{{suffix}}/metabat2",
        maxbin2 = f"{outdir}/binning/{{suffix}}/maxbin2",
        comebin = f"{outdir}/binning/{{suffix}}/comebin/comebin_res/comebin_res_bins",
        out_prefix = f"{outdir}/binning/{{suffix}}/binette"
    output:
        f"{outdir}/binning/{{suffix}}/binette/final_bins_quality_reports.tsv"
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

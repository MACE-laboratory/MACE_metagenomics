import pandas as pd

OUTDIR = config["outdir"]
RESOURCES = config["resources"]
BINNING_CONFIG = config["binning"]

threads_n = RESOURCES["threads"]
short_reads_folder = f"{OUTDIR}/trimmed_reads/short_reads"

groups = pd.read_csv(f"{OUTDIR}/assembly_groups.tsv", sep="\t", dtype=str)
suffixes = groups["assembly"].tolist()
assembly2samples = {
    row["assembly"]: [s.strip() for s in row["samples"].split(",") if s.strip()]
    for _, row in groups.iterrows()
}

rule all:
    input:
        expand(f"{OUTDIR}/assemblies/{{suffix}}/genes/eggnog/{{suffix}}_prokaryotic_genes.emapper.annotations", suffix=suffixes),
        expand(f"{OUTDIR}/assemblies/{{suffix}}/coverm/coverm_multi_sample.tsv", suffix=suffixes)

rule eggnog_mapper:
    input:
        assembly = f"{OUTDIR}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"
    output:
        ann = f"{OUTDIR}/assemblies/{{suffix}}/genes/eggnog/{{suffix}}_prokaryotic_genes.emapper.annotations"
    conda:
        "envs/eggnog_blast.yml"
    threads: threads_n
    shell:
        """
        mkdir -p {OUTDIR}/assemblies/{wildcards.suffix}/genes/eggnog
        emapper.py --cpu {threads} --dbmem --temp_dir /data/tmp_eggnog -m mmseqs --itype metagenome --data_dir /data/databases/eggnog-mapper2.1.2-2025 \
            -i {input.assembly} --resume \
            -o {OUTDIR}/assemblies/{wildcards.suffix}/genes/eggnog/{wildcards.suffix}_prokaryotic_genes
        """

def get_coverm_inputs(wildcards):
    sample_list = assembly2samples[wildcards.suffix]
    fq1s = [f"{short_reads_folder}/{sample}_R1_001_val_1.fq.gz" for sample in sample_list]
    fq2s = [f"{short_reads_folder}/{sample}_R2_001_val_2.fq.gz" for sample in sample_list]
    assembly = f"{OUTDIR}/assemblies/{wildcards.suffix}/prokaryotic_contigs.renamed.filtered.fa"
    return {"fq1s": fq1s, "fq2s": fq2s, "assembly": assembly}

rule coverm_multi:
    input:
        unpack(get_coverm_inputs)
    params:
        fq1=lambda wc, input: " ".join(input.fq1s),
        fq2=lambda wc, input: " ".join(input.fq2s)
    output:
        f"{OUTDIR}/assemblies/{{suffix}}/coverm/coverm_multi_sample.tsv"
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        """
        mkdir -p {OUTDIR}/assemblies/{wildcards.suffix}/coverm
        coverm contig \
            -t {threads} \
            -r {input.assembly} \
            -1 {params.fq1} \
            -2 {params.fq2} \
            -m trimmed_mean \
            --output-file {output}
        """

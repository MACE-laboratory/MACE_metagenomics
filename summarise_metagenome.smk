import pandas as pd

outdir = config["outdir"]
threads_n = config["threads"]
short_reads_folder = config["short_reads_folder"]

# Parse assembly_groups.tsv for suffix/sample mapping
groups = pd.read_csv(f"{outdir}/assembly_groups.tsv", sep="\t", dtype=str)
suffixes = groups["assembly"].tolist()
assembly2samples = {
    row["assembly"]: [s.strip() for s in row["samples"].split(",") if s.strip()]
    for _, row in groups.iterrows()
}

rule all:
    input:
        expand(f"{outdir}/assemblies/{{suffix}}/genes/eggnog/{{suffix}}_prokaryotic_genes.emapper.annotations", suffix=suffixes),
        expand(f"{outdir}/assemblies/{{suffix}}/coverm/coverm_multi_sample.tsv", suffix=suffixes)

rule eggnog_mapper:
    input:
        assembly = f"{outdir}/assemblies/{{suffix}}/prokaryotic_contigs.renamed.filtered.fa"
    output:
        ann = f"{outdir}/assemblies/{{suffix}}/genes/eggnog/{{suffix}}_prokaryotic_genes.emapper.annotations"
    conda:
        "eggnog_blast"
    threads: threads_n
    shell:
        """
        mkdir -p {outdir}/assemblies/{wildcards.suffix}/genes/eggnog
        emapper.py --cpu {threads} --dbmem --temp_dir /data/tmp_eggnog -m mmseqs --itype metagenome --data_dir /data/databases/eggnog-mapper2.1.2-2025 \
            -i {input.assembly} --resume \
            -o {outdir}/assemblies/{wildcards.suffix}/genes/eggnog/{wildcards.suffix}_prokaryotic_genes
        """

def get_coverm_inputs(wildcards):
    sample_list = assembly2samples[wildcards.suffix]
    fq1s = [f"{short_reads_folder}/{sample}_R1_001_val_1.fq.gz" for sample in sample_list]
    fq2s = [f"{short_reads_folder}/{sample}_R2_001_val_2.fq.gz" for sample in sample_list]
    assembly = f"{outdir}/assemblies/{wildcards.suffix}/prokaryotic_contigs.renamed.filtered.fa"
    return {"fq1s": fq1s, "fq2s": fq2s, "assembly": assembly}

rule coverm_multi:
    input:
        unpack(get_coverm_inputs)
    params:
        fq1=lambda wc, input: " ".join(input.fq1s),
        fq2=lambda wc, input: " ".join(input.fq2s)
    output:
        f"{outdir}/assemblies/{{suffix}}/coverm/coverm_multi_sample.tsv"
    conda:
        "envs/coverm.yml"
    threads: threads_n
    shell:
        """
        mkdir -p {outdir}/assemblies/{wildcards.suffix}/coverm
        coverm contig \
            -t {threads} \
            -r {input.assembly} \
            -1 {params.fq1} \
            -2 {params.fq2} \
            -m trimmed_mean \
            --output-file {output}
        """


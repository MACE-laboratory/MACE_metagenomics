import os
import datetime
import pandas as pd

################################################################################
# CONFIG (passed via --configfile in the snakemake call)
################################################################################

OUTDIR = config["outdir"]
CPUS = config["cpus"]
ASSEMBLIES_DIR = os.path.join(OUTDIR, "assemblies")
GROUPS_TSV = os.path.join(OUTDIR, "assembly_groups.tsv")

################################################################################
# Read assembly groups at parse time
################################################################################

ASSEMBLY_GROUPS = pd.read_csv(GROUPS_TSV, sep="\t")
ASSEMBLIES = ASSEMBLY_GROUPS["assembly"].tolist()
CONTIGS_DICT = dict(zip(ASSEMBLY_GROUPS["assembly"], ASSEMBLY_GROUPS["assembly_path"]))

################################################################################
# RULES
################################################################################

rule all:
    input:
        expand(os.path.join(ASSEMBLIES_DIR, "{assembly}/whokaryote/prokaryote_contig_headers.txt"), assembly=ASSEMBLIES),
        os.path.join(OUTDIR, "tools_version_log.txt")

rule run_whokaryote:
    input:
        contigs = lambda wildcards: CONTIGS_DICT[wildcards.assembly]
    output:
        outdir = directory(os.path.join(ASSEMBLIES_DIR, "{assembly}/whokaryote")),
        predictions = os.path.join(ASSEMBLIES_DIR, "{assembly}/whokaryote/prokaryote_contig_headers.txt")
    params:
        minsize = 1000
    threads:
        CPUS
    conda:
        "envs/whokaryote.yml"
    log:
        os.path.join(OUTDIR, "logs/whokaryote_{assembly}.log")
    shell:
        "whokaryote.py --contigs {input.contigs} --outdir {output.outdir} --minsize {params.minsize} --threads {threads}"

rule tools_version_log:
    input:
        expand(os.path.join(ASSEMBLIES_DIR, "{assembly}/whokaryote/prokaryote_contig_headers.txt"), assembly=ASSEMBLIES)
    output:
        os.path.join(OUTDIR, "tools_version_log.txt")
    run:
        with open(output[0], 'w') as log:
            log.write(f"Tools version log - {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            log.write("=================================\n\n")
            log.write("## Commands run:\n")
            log.write("=================================\n")
            for assembly in ASSEMBLIES:
                contigs = CONTIGS_DICT[assembly]
                outdir = os.path.join(ASSEMBLIES_DIR, assembly, "whokaryote")
                cmd = (f"whokaryote.py --contigs {contigs} --outdir {outdir} "
                       f"--minsize 1000 --threads {CPUS}")
                log.write(cmd + "\n")

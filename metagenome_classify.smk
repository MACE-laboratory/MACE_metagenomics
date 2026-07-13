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
        expand(os.path.join(ASSEMBLIES_DIR, "{assembly}/tiara_assignments.txt"), assembly=ASSEMBLIES),
        expand(os.path.join(ASSEMBLIES_DIR, "{assembly}/contigs_classification.done"), assembly=ASSEMBLIES),
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


rule run_tiara:
    input:
        contigs = lambda wildcards: CONTIGS_DICT[wildcards.assembly]
    output:
        os.path.join(ASSEMBLIES_DIR, "{assembly}/tiara_assignments.txt")
    params:
        minsize = 1000
    threads: 
        CPUS
    conda:
        "envs/tiara.yml"
    log: 
        os.path.join(OUTDIR, "logs/tiara_{assembly}.log")
    shell:
        "tiara -i {input.contigs} -m {params.minsize} -o {output} -t {threads}"

rule split_tiara_assignment_fasta:
    input:
        assignment = os.path.join(ASSEMBLIES_DIR, "{assembly}/tiara_assignments.txt"),
        contigs = lambda wildcards: CONTIGS_DICT[wildcards.assembly]
    output:
        prok = os.path.join(ASSEMBLIES_DIR, "{assembly}/prokaryotic_contigs.fasta"),
        euk = os.path.join(ASSEMBLIES_DIR, "{assembly}/eukaryotic_contigs.fasta"),
        other = os.path.join(ASSEMBLIES_DIR, "{assembly}/other_contigs.fasta"),
        done = os.path.join(ASSEMBLIES_DIR, "{assembly}/contigs_classification.done")
    params:
        assembly = "{assembly}"
    run:
        # Helper: Parse FASTA to a dict of {id: (header_line, sequence)}
        def parse_fasta(filename):
            seqs = {}
            current_id, current_header, seq_chunks = None, None, []
            with open(filename) as fh:
                for line in fh:
                    if line.startswith('>'):
                        if current_id is not None:
                            seqs[current_id] = (current_header, ''.join(seq_chunks))
                        current_header = line.strip()
                        current_id = line[1:].strip().split()[0]
                        seq_chunks = []
                    else:
                        seq_chunks.append(line.rstrip())
                # Don't forget last entry
                if current_id is not None:
                    seqs[current_id] = (current_header, ''.join(seq_chunks))
            return seqs

        # Read assignments
        id2class = {}
        with open(input.assignment) as f:
            header = f.readline()
            for line in f:
                if not line.strip():
                    continue
                seqid_col = line.strip().split('\t')[0]
                seqid = seqid_col.split()[0]
                klass = line.strip().split('\t')[1].lower()
                id2class[seqid] = klass

        prok_ids = {seqid for seqid, klass in id2class.items() if klass in ("bacteria", "archaea", "prokarya")}
        euk_ids = {seqid for seqid, klass in id2class.items() if klass == "eukarya"}
        other_ids = {seqid for seqid, klass in id2class.items() if klass not in ("bacteria", "archaea", "prokarya", "eukarya")}

        contigs = parse_fasta(input.contigs)

        # Helper: write out selected IDs
        def write_fasta(ids, fname):
            with open(fname, "w") as out:
                for seqid in ids:
                    if seqid in contigs:
                        header, seq = contigs[seqid]
                        # Prepend assembly name to ID in header
                        fields = header.split(maxsplit=1)
                        # Replace only first word after '>'
                        fields[0] = f">{params.assembly}_{seqid}"
                        out.write(' '.join(fields) + "\n")
                        # Write sequence in 60-char lines (as is customary)
                        for i in range(0, len(seq), 60):
                            out.write(seq[i:i+60] + "\n")

        write_fasta(prok_ids, output.prok)
        write_fasta(euk_ids, output.euk)
        write_fasta(other_ids, output.other)
        with open(output.done, "w") as f:
            f.write("done\n")

rule tools_version_log:
    input:
        expand(os.path.join(ASSEMBLIES_DIR, "{assembly}/tiara_assignments.txt"), assembly=ASSEMBLIES)
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
                cmd = (f"tiara -i {contigs} -o output -t {threads}")
                log.write(cmd + "\n")


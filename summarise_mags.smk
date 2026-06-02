import os
import glob

configfile: "config.yaml"

OUTDIR = config["outdir"]
SUFFIXES = config["suffixes"]

# ── Discover source bins from configured suffixes ────────────────────────────
#FA_FILES = []
#for suffix in SUFFIXES:
#    pattern = os.path.join(OUTDIR, "binning", suffix, "binette", "final_bins", "*.fa")
#    FA_FILES.extend(glob.glob(pattern))

#FA_FILES = sorted(FA_FILES)

#if not FA_FILES:
#    searched = [
#        os.path.join(OUTDIR, "binning", suffix, "binette", "final_bins", "*.fa")
#        for suffix in SUFFIXES
#    ]
#    raise ValueError(
#        "No .fa files found in any configured suffix directory. Searched:\n  - "
#        + "\n  - ".join(searched)
#    )

# ── Discover trimmed reads at parse time ──────────────────────────────────────
FWD_READS = sorted(
    glob.glob(os.path.join(OUTDIR, "trimmed_reads", "short_reads", "*_R1_001_val_1.fq.gz"))
)
REV_READS = sorted(
    glob.glob(os.path.join(OUTDIR, "trimmed_reads", "short_reads", "*_R2_001_val_2.fq.gz"))
)

if not FWD_READS or not REV_READS:
    raise ValueError(
        f"No trimmed read pairs found under {OUTDIR}/trimmed_reads/short_reads/."
    )

if len(FWD_READS) != len(REV_READS):
    raise ValueError(
        f"Mismatched number of forward ({len(FWD_READS)}) and reverse ({len(REV_READS)}) read files."
    )

# ── Path helpers ──────────────────────────────────────────────────────────────
MAGS_DIR                 = os.path.join(OUTDIR, "MAGs")

# pre-dereplication
ALL_BINS_DIR             = os.path.join(MAGS_DIR, "all_bins")
ALL_BINS_CHECKM2_DIR     = os.path.join(ALL_BINS_DIR, "checkm2")

# post-dereplication / final MAG dataset
GALAH_DIR                = os.path.join(MAGS_DIR, "galah")
MAGS_FASTAS_DIR          = os.path.join(MAGS_DIR, "mags_fastas")
MAGS_CHECKM2_DIR         = os.path.join(MAGS_DIR, "checkm2")
GTDBTK_DIR               = os.path.join(MAGS_DIR, "gtdb")
COVERM_DIR               = os.path.join(MAGS_DIR, "coverm")
#EMAPPER_DIR              = os.path.join(MAGS_DIR, "emapper")
BAKTA_DIR                = os.path.join(MAGS_DIR, "bakta")
ALL_MAGS_DIR             = os.path.join(MAGS_DIR, "all_mags")

# helper outputs
REP_LIST                 = os.path.join(GALAH_DIR, "rep_list")
GALAH_CLUSTER_TSV        = os.path.join(GALAH_DIR, "dereplicated_genomes.tsv")
GALAH_DONE               = os.path.join(GALAH_DIR, ".galah_done")
MAGS_FASTAS_DONE         = os.path.join(MAGS_FASTAS_DIR, ".mags_fastas_done")
MAGS_QUALITY_REPORT      = os.path.join(MAGS_CHECKM2_DIR, "mags_quality_report.tsv")
GTDBTK_DONE              = os.path.join(GTDBTK_DIR, "gtdbtk.done")
#EMAPPER_DONE             = os.path.join(EMAPPER_DIR, ".emapper_done")
BAKTA_DONE             = os.path.join(BAKTA_DIR, ".bakta_done")
COVERM_TSV               = os.path.join(COVERM_DIR, "mags_abundances.tsv")

# ── Helper functions ──────────────────────────────────────────────────────────
def get_mag_bases():
    fastas = sorted(glob.glob(os.path.join(MAGS_FASTAS_DIR, "*.fa")))
    return [os.path.splitext(os.path.basename(f))[0] for f in fastas]


def get_protein_outputs(wildcards):
    return expand(
        os.path.join(ALL_MAGS_DIR, "{base}.faa"),
        base=get_mag_bases()
    )


def get_emapper_outputs(wildcards):
    return expand(
        os.path.join(EMAPPER_DIR, "results", "{base}.emapper.annotations"),
        base=get_mag_bases()
    )

def get_bakta_outputs(wildcards):
    return expand(
        os.path.join(BAKTA_DIR, "results", "{base}", "{base}.tsv"),
        base=get_mag_bases()
    )

# ── Rules ─────────────────────────────────────────────────────────────────────

rule all:
    input:
        os.path.join(ALL_BINS_DIR, ".all_bins_copied"),
        os.path.join(ALL_BINS_CHECKM2_DIR, "quality_report.tsv"),
        GALAH_CLUSTER_TSV,
        GALAH_DONE,
        MAGS_FASTAS_DONE,
        MAGS_QUALITY_REPORT,
        GTDBTK_DONE,
        COVERM_TSV,
        BAKTA_DONE,        
#EMAPPER_DONE,


# ── 1. Copy all bins to MAGs/all_bins/ ───────────────────────────────────────

rule copy_all_bins:
    output:
        sentinel=os.path.join(ALL_BINS_DIR, ".all_bins_copied"),
    params:
        outdir=ALL_BINS_DIR,
        outdir_root=OUTDIR,
        suffixes=SUFFIXES,
    run:
        import os
        import glob
        import shutil

        os.makedirs(params.outdir, exist_ok=True)

        fa_files = []
        searched_patterns = []

        for suffix in params.suffixes:
            pattern = os.path.join(
                params.outdir_root,
                "binning",
                suffix,
                "binette",
                "final_bins",
                "*.fa"
            )
            searched_patterns.append(pattern)

            for src in sorted(glob.glob(pattern)):
                fa_files.append((suffix, src))

        if not fa_files:
            raise ValueError(
                "No .fa files found for copy_all_bins. Searched:\n  - "
                + "\n  - ".join(searched_patterns)
            )

        # remove stale copied fasta files from previous runs
        for old_fa in glob.glob(os.path.join(params.outdir, "*.fa")):
            os.remove(old_fa)

        copied = []
        seen_dest_names = set()

        for suffix, src in fa_files:
            bin_name = os.path.basename(src)
            dest_name = f"{suffix}_{bin_name}"
            dest = os.path.join(params.outdir, dest_name)

            if dest_name in seen_dest_names:
                raise ValueError(
                    f"Duplicate destination filename after renaming: {dest_name}"
                )

            shutil.copy2(src, dest)
            seen_dest_names.add(dest_name)
            copied.append(dest)

        with open(output.sentinel, "w") as fh:
            fh.write("\n".join(copied) + "\n")

# ── 2. CheckM2 on all bins before dereplication ───────────────────────────────

rule checkm2_all_bins:
    input:
        sentinel=os.path.join(ALL_BINS_DIR, ".all_bins_copied"),
    output:
        tsv=os.path.join(ALL_BINS_CHECKM2_DIR, "quality_report.tsv"),
    params:
        input_dir=ALL_BINS_DIR,
        outdir=ALL_BINS_CHECKM2_DIR,
        extension="fa",
    threads: 32
    conda:
        "checkm2"
    log:
        os.path.join(OUTDIR, "logs", "checkm2", "checkm2_all_bins.log")
    shell:
        """
        mkdir -p {params.outdir}

        checkm2 predict \
            --threads {threads} \
            --input {params.input_dir} \
            -x {params.extension} \
            --output-directory {params.outdir} \
        2>&1 | tee {log}
        """


# ── 3. Galah dereplication ────────────────────────────────────────────────────

rule galah_dereplicate:
    input:
        bins_done=os.path.join(ALL_BINS_DIR, ".all_bins_copied"),
        checkm2_tsv=os.path.join(ALL_BINS_CHECKM2_DIR, "quality_report.tsv"),
    output:
        tsv=GALAH_CLUSTER_TSV,
        rep_list=REP_LIST,
        sentinel=GALAH_DONE,
    params:
        genome_dir=ALL_BINS_DIR,
        outdir=GALAH_DIR,
        rep_dir=os.path.join(GALAH_DIR, "cluster_rep"),
        fasta_copy=os.path.join(GALAH_DIR, "fasta_copy"),
        extension="fa",
    threads: 32
    conda:
        "envs/galah.yml"
    log:
        os.path.join(OUTDIR, "logs", "galah", "galah_dereplicate.log")
    shell:
        """
        mkdir -p {params.outdir} {params.rep_dir} {params.fasta_copy}

        galah cluster \
            --genome-fasta-directory {params.genome_dir} \
            --genome-fasta-extension {params.extension} \
            --checkm2-quality-report {input.checkm2_tsv} \
            --ani 99 \
            --min-completeness 40 \
            --max-contamination 20 \
            -t {threads} \
            --output-cluster-definition {output.tsv} \
            --output-representative-fasta-directory {params.rep_dir} \
            --output-representative-fasta-directory-copy {params.fasta_copy} \
            --output-representative-list {output.rep_list} \
        2>&1 | tee {log}

        touch {output.sentinel}
        """


# ── 4. Copy dereplicated fasta set into MAGs/mags_fastas/ ────────────────────

rule collect_mags_fastas:
    input:
        sentinel=GALAH_DONE,
    output:
        sentinel=MAGS_FASTAS_DONE,
    params:
        src_dir=os.path.join(GALAH_DIR, "cluster_rep"),
        outdir=MAGS_FASTAS_DIR,
    shell:
        """
        mkdir -p {params.outdir}
        find {params.outdir} -maxdepth 1 -type f -name "*.fa" -delete

        for f in {params.src_dir}/*.fa; do
            [ -e "$f" ] || continue
            cp "$f" {params.outdir}/
        done

        touch {output.sentinel}
        """


# ── 5. Subset CheckM2 report to retained MAGs only ───────────────────────────

rule subset_checkm2_for_mags:
    input:
        mags_done=MAGS_FASTAS_DONE,
        checkm2_tsv=os.path.join(ALL_BINS_CHECKM2_DIR, "quality_report.tsv"),
    output:
        tsv=MAGS_QUALITY_REPORT,
    params:
        mags_dir=MAGS_FASTAS_DIR,
        outdir=MAGS_CHECKM2_DIR,
    log:
        os.path.join(OUTDIR, "logs", "checkm2", "subset_mags_quality_report.log")
    run:
        import os
        import pandas as pd

        # Get MAG names from fasta files in MAGS_FASTA_DIR, stripping ".fa"
        mag_names = {
            os.path.splitext(f)[0]
            for f in os.listdir(params.mags_dir)
            if f.endswith(".fa")
        }

        # Load CheckM2 report
        df = pd.read_csv(input.checkm2_tsv, sep="\t")

        # Subset rows where Name matches MAG names
        df_subset = df[df["Name"].isin(mag_names)]

        # Write subsetted TSV
        df_subset.to_csv(output.tsv, sep="\t", index=False) 
        
        

# ── 6. GTDB-Tk on dereplicated MAGs ───────────────────────────────────────────

rule gtdbtk:
    input:
        mags_done=MAGS_FASTAS_DONE,
    output:
        done=GTDBTK_DONE,
    params:
        mag_dir=MAGS_FASTAS_DIR,
        outdir=GTDBTK_DIR,
        extension="fa",
        gtdb_data="/data/databases/gtdb/release232/",
    threads: 32
    conda:
        "gtdbtk-2.7.2"
    log:
        os.path.join(OUTDIR, "logs", "gtdbtk", "gtdbtk.log")
    shell:
        """
        mkdir -p {params.outdir}

        export GTDBTK_DATA_PATH="{params.gtdb_data}"

        gtdbtk classify_wf \
            --genome_dir {params.mag_dir} \
            -x {params.extension} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            # that was for previous version --skip_ani_screen \
        2>&1 | tee {log}

        touch {output.done}
        """


# ── 7. CoverM on dereplicated MAGs ────────────────────────────────────────────

rule coverm:
    input:
        mags_done=MAGS_FASTAS_DONE,
        fwd=FWD_READS,
        rev=REV_READS,
    output:
        tsv=COVERM_TSV,
    params:
        mag_dir=MAGS_FASTAS_DIR,
        extension="fa",
        fwd_str=lambda wildcards, input: " ".join(input.fwd),
        rev_str=lambda wildcards, input: " ".join(input.rev),
    threads: 64
    conda:
        "envs/coverm.yml"
    log:
        os.path.join(OUTDIR, "logs", "coverm", "coverm.log")
    shell:
        """
        mkdir -p {COVERM_DIR}

        coverm genome \
            --genome-fasta-directory {params.mag_dir} \
            --genome-fasta-extension {params.extension} \
            -1 {params.fwd_str} \
            -2 {params.rev_str} \
            -p bwa-mem2 \
            -m trimmed_mean \
            -o {output.tsv} \
            -t {threads} \
        2>&1 | tee {log}
        """


# ── 8. Prodigal on each MAG fasta ─────────────────────────────────────────────

#rule prodigal:
#    input:
#        mag=os.path.join(MAGS_FASTAS_DIR, "{base}.fa"),
#        mags_done=MAGS_FASTAS_DONE,
#    output:
#        faa=os.path.join(ALL_MAGS_DIR, "{base}.faa"),
#    threads: 1
#    conda:
#        "eggnog_blast"
#    log:
#        os.path.join(OUTDIR, "logs", "emapper", "prodigal_{base}.log")
#    shell:
#        """
#        mkdir -p {ALL_MAGS_DIR}
#
#        prodigal \
#            -i {input.mag} \
#            -a {output.faa} \
#            -p meta \
#        2>&1 | tee {log}
#        """


# ── 9. eggNOG-mapper on each MAG protein set ─────────────────────────────────

#rule emapper:
#    input:
#        faa=os.path.join(ALL_MAGS_DIR, "{base}.faa"),
#    output:
#        annot=os.path.join(EMAPPER_DIR, "results", "{base}.emapper.annotations"),
#    params:
#        outdir=os.path.join(EMAPPER_DIR, "results"),
#        prefix=lambda wildcards: wildcards.base,
#    threads: 32
#    conda:
#        "eggnog_blast"
#    log:
#        os.path.join(OUTDIR, "logs", "emapper", "emapper_{base}.log")
#    shell:
#        """
#        mkdir -p {params.outdir}
#
#        emapper.py \
#            -i {input.faa} \
#            -o {params.prefix} \
#            --output_dir {params.outdir} \
#            --cpu {threads} \
#            -m mmseqs \
#            --itype proteins \
#        2>&1 | tee {log}
#        """

# bakta too

rule bakta:
    input:
        fasta=os.path.join(MAGS_FASTAS_DIR, "{base}.fa"),
    output:
        annot=os.path.join(BAKTA_DIR, "results", "{base}", "{base}.tsv"),
    params:
        outdir=lambda wildcards: os.path.join(BAKTA_DIR, wildcards.base),
        prefix=lambda wildcards: wildcards.base,
        gz=lambda wildcards: os.path.join(BAKTA_DIR, "gzipped_fastas", wildcards.base, f"{wildcards.base}.fa.gz"),
    threads: 8
    conda:
        "envs/bakta.yml"
    log:
        os.path.join(OUTDIR, "logs", "bakta", "bakta_{base}.log")
    shell:
        """
        mkdir -p $(dirname {params.gz})

        # gzip input fasta for Bakta
        gzip -c {input.fasta} > {params.gz}

        bakta \
            --db /data/databases/bakta/bakta_full_01062026/db \
            --output {params.outdir} \
            --prefix {params.prefix} \
            --threads {threads} \
            {params.gz} \
        2>&1 | tee {log}
        """

# ── 10. Final emapper sentinel ────────────────────────────────────────────────

#rule emapper_done:
#    input:
#        get_emapper_outputs,
#    output:
#        sentinel=EMAPPER_DONE,
#    shell:
#        """
#        touch {output.sentinel}
#        """

#bakta
rule bakta_done:
    input:
        get_bakta_outputs,
    output:
        sentinel=BAKTA_DONE,
    shell:
        """
        touch {output.sentinel}
        """

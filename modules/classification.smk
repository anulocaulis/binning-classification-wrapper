# modules/classification.smk
# Snakemake module for taxonomic classification of reads and MAGs
# Tools: Kraken2 (read classification), GTDB-Tk (bin classification)

import os


def resolve_kraken2_db():
    configured = config.get("kraken2_db", "")
    candidates = [configured]

    if configured.endswith("/kraken2_db"):
        candidates.append(configured[:-len("/kraken2_db")] + "/k2_standard")

    for candidate in config.get("kraken2_db_fallbacks", []):
        candidates.append(candidate)

    for db in candidates:
        if db and all(os.path.exists(os.path.join(db, f)) for f in ("hash.k2d", "taxo.k2d", "opts.k2d")):
            return db

    return configured

# ---
# Step 1: Kraken2 Read Classification
# ---
rule kraken2_classify:
    """
    Classifies metagenomic reads using Kraken2 against the configured database.
    Produces a classification report and per-read output.
    """
    input:
        reads = lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample)
    output:
        report = f"{OUTPUT_DIR}/{{sample}}/classification/kraken2/report.txt",
        output_file = f"{OUTPUT_DIR}/{{sample}}/classification/kraken2/output.kraken"
    params:
        outdir = f"{OUTPUT_DIR}/{{sample}}/classification/kraken2",
        container = CLASSIFICATION_CONTAINER,
        db = resolve_kraken2_db()
    threads: config["threads"]
    log: "logs/classification_kraken2_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        test -f {params.db}/hash.k2d
        test -f {params.db}/taxo.k2d
        test -f {params.db}/opts.k2d
        singularity exec {params.container} /opt/conda/envs/base_tools/bin/kraken2 \
            --db {params.db} \
            --threads {threads} \
            --report {output.report} \
            --output {output.output_file} \
            {input.reads} 2>> {log}
        """

# ---
# Step 2: GTDB-Tk Bin Classification
# ---
rule gtdbtk_classify:
    """
    Classifies refined MAG bins using GTDB-Tk to assign taxonomy based on the
    GTDB reference database.
    """
    input:
        bins_dir = f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
    output:
        summary = f"{OUTPUT_DIR}/{{sample}}/classification/gtdbtk/gtdbtk.bac120.summary.tsv"
    params:
        outdir = f"{OUTPUT_DIR}/{{sample}}/classification/gtdbtk",
        container = CLASSIFICATION_CONTAINER,
        db = config["gtdbtk_db"]
    threads: config["threads"]
    log: "logs/classification_gtdbtk_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} gtdbtk classify_wf \
            --genome_dir {input.bins_dir} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --extension .fa \
            --skip_ani_screen 2>> {log}
        """

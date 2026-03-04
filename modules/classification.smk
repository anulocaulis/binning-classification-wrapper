# modules/classification.smk
# Snakemake module for taxonomic classification of reads and MAGs
# Tools: Kraken2 (read classification), GTDB-Tk (bin classification)

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
        report = "{output_dir}/{sample}/classification/kraken2/report.txt",
        output_file = "{output_dir}/{sample}/classification/kraken2/output.kraken"
    params:
        outdir = "{output_dir}/{sample}/classification/kraken2",
        container = CLASSIFICATION_CONTAINER,
        db = config["kraken2_db"]
    threads: config["threads"]
    log: "logs/classification_kraken2_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} kraken2 \
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
        bin_refinement = "{output_dir}/{sample}/bin_refinement"
    output:
        summary = "{output_dir}/{sample}/classification/gtdbtk/gtdbtk.bac120.summary.tsv"
    params:
        outdir = "{output_dir}/{sample}/classification/gtdbtk",
        bins_dir = "{output_dir}/{sample}/bin_refinement/metawrap_50_10_bins",
        container = CLASSIFICATION_CONTAINER,
        db = config["gtdbtk_db"]
    threads: config["threads"]
    log: "logs/classification_gtdbtk_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --extension .fa \
            --skip_ani_screen 2>> {log}
        """

# ---
# Step 1: Kraken2 Taxonomic Classification
# ---
rule kraken2_classify:
    """
    Runs Kraken2 to assign taxonomic labels to reads using a k-mer based
    approach against the specified database. Produces a full classification
    output file and a summary report.
    """
    input:
        reads = config["input_reads"]["short_interleaved"]
    output:
        report = "{classification_output_dir}/{sample}/kraken2/report.txt",
        output = "{classification_output_dir}/{sample}/kraken2/output.txt"
    params:
        container = CLASSIFICATION_CONTAINER,
        db = config["kraken2_db"],
        outdir = "{classification_output_dir}/{sample}/kraken2"
    threads: config["threads"]
    log: "logs/classification_kraken2_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} kraken2 \
            --db {params.db} \
            --threads {threads} \
            --report {output.report} \
            --output {output.output} \
            {input.reads} 2>> {log}
        """


# ---
# Step 2: Bracken Abundance Re-estimation
# ---
rule bracken_abundance:
    """
    Runs Bracken to re-estimate species-level abundances from Kraken2
    classification report using a Bayesian approach.
    """
    input:
        report = "{classification_output_dir}/{sample}/kraken2/report.txt"
    output:
        bracken_report = "{classification_output_dir}/{sample}/bracken/{sample}_bracken.txt"
    params:
        container = CLASSIFICATION_CONTAINER,
        db = config["kraken2_db"],
        outdir = "{classification_output_dir}/{sample}/bracken"
    log: "logs/classification_bracken_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.bracken_report} 2>> {log}
        """


# ---
# Step 3: MetaPhlAn Taxonomic Profiling
# ---
rule metaphlan_profile:
    """
    Runs MetaPhlAn to profile the taxonomic composition of a metagenome using
    marker genes. Complements Kraken2/Bracken with a marker-gene based approach.
    """
    input:
        reads = config["input_reads"]["short_interleaved"]
    output:
        profile = "{classification_output_dir}/{sample}/metaphlan/{sample}_metaphlan_profile.txt"
    params:
        container = CLASSIFICATION_CONTAINER,
        outdir = "{classification_output_dir}/{sample}/metaphlan"
    threads: config["threads"]
    log: "logs/classification_metaphlan_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} metaphlan \
            {input.reads} \
            --input_type fastq \
            --nproc {threads} \
            -o {output.profile} 2>> {log}
        """

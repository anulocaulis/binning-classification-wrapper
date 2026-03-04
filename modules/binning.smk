# modules/binning.smk
# Snakemake module for metagenome-assembled genome (MAG) binning
# Tools: MetaWRAP (MetaBAT2, MaxBin2, CONCOCT), bin refinement, CheckM2

# ---
# Step 1: MetaWRAP Binning (MetaBAT2, MaxBin2, CONCOCT)
# ---
rule metawrap_binning:
    """
    Runs MetaWRAP binning module using MetaBAT2, MaxBin2, and CONCOCT to generate
    initial bin sets from the assembled contigs and mapped reads.
    """
    input:
        assembly = lambda wildcards: config["input_reads"]["assembly"].format(sample=wildcards.sample),
        reads = lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample)
    output:
        metabat2 = directory("{output_dir}/{sample}/binning/metabat2"),
        maxbin2 = directory("{output_dir}/{sample}/binning/maxbin2"),
        concoct = directory("{output_dir}/{sample}/binning/concoct")
    params:
        outdir = "{output_dir}/{sample}/binning",
        container = BINNING_CONTAINER,
        min_contig_len = 1000
    threads: config["threads"]
    log: "logs/binning_metawrap_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} metawrap binning \
            -o {params.outdir} \
            -t {threads} \
            -a {input.assembly} \
            --metabat2 --maxbin2 --concoct \
            -l {params.min_contig_len} \
            {input.reads} 2>> {log}
        """

# ---
# Step 2: MetaWRAP Bin Refinement
# ---
rule metawrap_bin_refinement:
    """
    Refines and consolidates bins from multiple binning algorithms using MetaWRAP's
    bin_refinement module. Produces a high-quality, non-redundant bin set.
    """
    input:
        metabat2 = "{output_dir}/{sample}/binning/metabat2",
        maxbin2 = "{output_dir}/{sample}/binning/maxbin2",
        concoct = "{output_dir}/{sample}/binning/concoct"
    output:
        refined_bins = directory("{output_dir}/{sample}/bin_refinement")
    params:
        outdir = "{output_dir}/{sample}/bin_refinement",
        refined_bins_dir = "{output_dir}/{sample}/bin_refinement/metawrap_50_10_bins",
        container = BINNING_CONTAINER,
        completeness = 50,
        contamination = 10
    threads: config["threads"]
    log: "logs/binning_refinement_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} metawrap bin_refinement \
            -o {params.outdir} \
            -t {threads} \
            -A {input.metabat2} \
            -B {input.maxbin2} \
            -C {input.concoct} \
            -c {params.completeness} \
            -x {params.contamination} 2>> {log}
        """

# ---
# Step 3: CheckM2 Bin Quality Assessment
# ---
rule checkm2_quality:
    """
    Assesses the completeness and contamination of refined bins using CheckM2.
    Produces a quality report for all bins.
    """
    input:
        bin_refinement = "{output_dir}/{sample}/bin_refinement"
    output:
        report = "{output_dir}/{sample}/checkm2/quality_report.tsv"
    params:
        outdir = "{output_dir}/{sample}/checkm2",
        bins_dir = "{output_dir}/{sample}/bin_refinement/metawrap_50_10_bins",
        container = BINNING_CONTAINER
    threads: config["threads"]
    log: "logs/binning_checkm2_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} checkm2 predict \
            --input {params.bins_dir} \
            --output-directory {params.outdir} \
            --threads {threads} \
            --extension .fa 2>> {log}
        """

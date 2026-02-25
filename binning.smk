# --- Helper Function for Assembly Input ---
def get_assembly_file(wildcards):
    """
    Returns the path to the assembled FASTA file for a given sample.
    Assemblies are expected to come from the assembly-eval-wrapper output.
    """
    return f"{config['assemblies']['dir']}/{wildcards.sample}/assembly.fasta"


# ---
# Step 1: Coverage Depth for MetaBAT2
# ---
rule metabat2_depth:
    """
    Computes per-contig coverage depth from a sorted BAM file using MetaBAT2's
    jgi_summarize_bam_contig_depths helper. This depth file is required for
    MetaBAT2 binning.
    """
    input:
        bam = "mapping/{sample}.sorted.bam"
    output:
        depth = "mapping/{sample}_depth.txt"
    params:
        container = METAWRAP_CONTAINER
    log: "logs/binning_metabat2_depth_{sample}.log"
    shell:
        "singularity exec {params.container} jgi_summarize_bam_contig_depths "
        "--outputDepth {output.depth} {input.bam} 2>> {log}"


# ---
# Step 2: MetaBAT2 Binning
# ---
rule metabat2_binning:
    """
    Runs MetaBAT2 to bin assembled contigs using coverage depth information.
    Produces a set of FASTA bin files in the output directory.
    """
    input:
        assembly = get_assembly_file,
        depth = "mapping/{sample}_depth.txt"
    output:
        done = touch("{output_dir}/{sample}/metabat2/.done")
    params:
        container = METAWRAP_CONTAINER,
        outdir = "{output_dir}/{sample}/metabat2"
    threads: config["threads"]
    log: "logs/binning_metabat2_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} metabat2 \
            -i {input.assembly} \
            -a {input.depth} \
            -o {params.outdir}/bin \
            -t {threads} 2>> {log}
        """


# ---
# Step 3: MaxBin2 Binning
# ---
rule maxbin2_binning:
    """
    Runs MaxBin2 to bin assembled contigs using read abundance information.
    Produces a set of FASTA bin files in the output directory.
    """
    input:
        assembly = get_assembly_file,
        reads = config["input_reads"]["short_interleaved"]
    output:
        done = touch("{output_dir}/{sample}/maxbin2/.done")
    params:
        container = METAWRAP_CONTAINER,
        outdir = "{output_dir}/{sample}/maxbin2"
    threads: config["threads"]
    log: "logs/binning_maxbin2_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} run_MaxBin.pl \
            -contig {input.assembly} \
            -reads {input.reads} \
            -out {params.outdir}/bin \
            -thread {threads} 2>> {log}
        """


# ---
# Step 4: Bin Quality Assessment with CheckM2
# ---
rule checkm2_evaluation:
    """
    Runs CheckM2 to evaluate the quality (completeness and contamination) of
    genome bins produced by a binning tool.
    """
    input:
        done = "{output_dir}/{sample}/{binner}/.done"
    output:
        report = "{output_dir}/{sample}/{binner}/checkm2/quality_report.tsv"
    params:
        container = METAWRAP_CONTAINER,
        bins_dir = "{output_dir}/{sample}/{binner}",
        outdir = "{output_dir}/{sample}/{binner}/checkm2"
    threads: config["threads"]
    log: "logs/binning_checkm2_{binner}_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        singularity exec {params.container} checkm2 predict \
            --input {params.bins_dir} \
            --output-directory {params.outdir} \
            --threads {threads} 2>> {log}
        """

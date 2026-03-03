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
        mkdir -p {params.outdir}/bins
        singularity exec {params.container} metabat2 \
            -i {input.assembly} \
            -a {input.depth} \
            -o {params.outdir}/bins/bin \
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
        mkdir -p {params.outdir}/bins
        singularity exec {params.container} run_MaxBin.pl \
            -contig {input.assembly} \
            -reads {input.reads} \
            -out {params.outdir}/bins/bin \
            -thread {threads} 2>> {log}
        """


# ---
# Step 3b: CONCOCT Binning
# ---
rule concoct_binning:
    """
    Runs CONCOCT to bin assembled contigs using a Gaussian mixture model on
    sequence composition and coverage. Requires a sorted BAM for coverage
    profiling.
    """
    input:
        assembly = get_assembly_file,
        bam = "mapping/{sample}.sorted.bam"
    output:
        done = touch("{output_dir}/{sample}/concoct/.done")
    params:
        container = METAWRAP_CONTAINER,
        outdir = "{output_dir}/{sample}/concoct"
    threads: config["threads"]
    log: "logs/binning_concoct_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}/cut_contigs {params.outdir}/bins

        # 1. Cut assembly into 10 kb chunks for coverage estimation
        singularity exec {params.container} cut_up_fasta.py {input.assembly} \
            -c 10000 -o 0 --merge_last \
            -b {params.outdir}/contigs_10K.bed \
            > {params.outdir}/contigs_10K.fa 2>> {log}

        # 2. Generate per-contig coverage table from BAM
        singularity exec {params.container} concoct_coverage_table.py \
            {params.outdir}/contigs_10K.bed {input.bam} \
            > {params.outdir}/coverage_table.tsv 2>> {log}

        # 3. Run CONCOCT clustering
        singularity exec {params.container} concoct \
            --composition_file {params.outdir}/contigs_10K.fa \
            --coverage_file {params.outdir}/coverage_table.tsv \
            -b {params.outdir}/concoct_output/ \
            -t {threads} 2>> {log}

        # 4. Merge cut-up clustering back to original contigs
        singularity exec {params.container} merge_cutup_clustering.py \
            {params.outdir}/concoct_output/clustering_gt1000.csv \
            > {params.outdir}/clustering_merged.csv 2>> {log}

        # 5. Extract individual bin FASTA files
        singularity exec {params.container} extract_fasta_bins.py \
            {input.assembly} {params.outdir}/clustering_merged.csv \
            --output_path {params.outdir}/bins/ 2>> {log}
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
        bins_dir = "{output_dir}/{sample}/{binner}/bins",
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


# ---
# Step 5: Bin Taxonomic Classification with GTDB-Tk
# ---
rule gtdbtk_classify:
    """
    Runs GTDB-Tk to assign standardised taxonomy to genome bins using the
    GTDB reference database. Requires CheckM2 to have completed first so that
    only high-quality bins are classified.
    """
    input:
        done = "{output_dir}/{sample}/{binner}/.done"
    output:
        summary = "{output_dir}/{sample}/{binner}/gtdbtk/gtdbtk.bac120.summary.tsv"
    params:
        container = METAWRAP_CONTAINER,
        bins_dir = "{output_dir}/{sample}/{binner}/bins",
        outdir = "{output_dir}/{sample}/{binner}/gtdbtk",
        gtdbtk_db = config["gtdbtk_db"]
    threads: config["threads"]
    log: "logs/binning_gtdbtk_{binner}_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        GTDBTK_DATA_PATH={params.gtdbtk_db} \
        singularity exec {params.container} gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir {params.outdir} \
            --extension fa \
            --cpus {threads} 2>> {log}
        """

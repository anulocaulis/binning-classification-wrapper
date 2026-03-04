

configfile: "config.yaml"

ASSEMBLY_CONTAINER = "containers/assembler.sif"
QC_CONTAINER = "containers/qc_tools_miniconda.sif"
FLYE_ASSEMBLY_CONTAINER = "containers/flye_assembler.sif"
NANOPLOT_CONTAINER = "containers/nano-tools4.sif"
POLISHER_CONTAINER = "containers/nextpolish.sif"

NanoPlot_tools = ["flye", "metaspades", "nextpolish"]

rule all:
    input:
        expand("{output_dir}/{sample}/assembly.flye/assembly.fasta", output_dir=config["output_dir"], sample=config["long_read_samples"]),
        expand("{output_dir}/{sample}/assembly.metamdbg/contigs.fasta", output_dir=config["output_dir"], sample=config["long_read_samples"]),
        expand("data/{sample}_long_reads_filtered.fastq.gz", sample=config["long_read_samples"]),
        # Disabled - IDBA-UD cannot handle gzipped FASTQ input, requires FASTA conversion
        # expand("{output_dir}/{sample}/assembly.idbaud/assembly.fasta", output_dir=config["output_dir"], sample=config["short_read_samples"]),
        expand("{output_dir}/{sample}/assembly.metaspades/assembly.fasta", output_dir=config["output_dir"], sample=config["short_read_samples"]),
        expand("{output_dir}/{sample}/assembly.megahit/final.contigs.fa", output_dir=config["output_dir"], sample=config["short_read_samples"]),
        expand("{output_dir}/{sample}/assembly.metaspades_hybrid/assembly.fasta", output_dir=config["output_dir"], sample=config["hybrid_samples"]),
        expand("{output_dir}/{sample}/assembly.nextpolish/assembly.fasta", output_dir=config["output_dir"], sample=config["hybrid_samples"]),
        expand("qc/{tool}_{sample}/quast/report.html", tool=["flye","metamdbg"], sample=config["long_read_samples"]),
        expand("qc/{tool}_{sample}/quast/report.html", tool=["metaspades","megahit"], sample=config["short_read_samples"]),
        expand("qc/{tool}_{sample}/quast/report.html", tool=["metaspades_hybrid","nextpolish"], sample=config["hybrid_samples"]),
        expand("qc/{tool}_{sample}/assembly_stats.txt", tool=["flye","metamdbg"], sample=config["long_read_samples"]),
        expand("qc/{tool}_{sample}/assembly_stats.txt", tool=["metaspades","megahit"], sample=config["short_read_samples"]),
        expand("qc/{tool}_{sample}/assembly_stats.txt", tool=["metaspades_hybrid","nextpolish"], sample=config["hybrid_samples"]),
        # Temporarily disabled bbmap_qc due to timeouts
        # expand("qc/{tool}_{sample}/bbmap/stats_summary.txt", tool=["flye","metamdbg"], sample=config["long_read_samples"]),
        # expand("qc/{tool}_{sample}/bbmap/stats_summary.txt", tool=["idbaud","metaspades","megahit"], sample=config["short_read_samples"]),
        # expand("qc/{tool}_{sample}/bbmap/stats_summary.txt", tool=["metaspades_hybrid","nextpolish"], sample=config["hybrid_samples"]),
        # Disabled - these FastQC paths don't match any rules in qc.smk
        # expand("qc/{tool}_{sample}/fastqc/{sample}_interleaved_fastqc.zip", tool=["flye","metamdbg"], sample=config["long_read_samples"]),
        # expand("qc/{tool}_{sample}/fastqc/{sample}_interleaved_fastqc.zip", tool=["idbaud","metaspades","megahit"], sample=config["short_read_samples"]),
        # expand("qc/{tool}_{sample}/fastqc/{sample}_interleaved_fastqc.zip", tool=["metaspades_hybrid","nextpolish"], sample=config["hybrid_samples"]),
        expand("qc/{tool}_{sample}/nanoplot/NanoPlot-report.html", tool=NanoPlot_tools, sample=config["long_read_samples"]),
        expand("qc/{tool}_{sample}/nanoqc/nanoQC.html", tool=NanoPlot_tools, sample=config["long_read_samples"])

include: "modules/read_prep.smk"
include: "modules/assembly.smk"
include: "modules/qc.smk"

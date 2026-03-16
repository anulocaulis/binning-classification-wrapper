# modules/preassembly_qc.smk
# Pre-assembly quality control for Illumina and ONT inputs

rule fastqc_short_reads:
    """
    Runs FastQC on interleaved short-read FASTQ input.
    """
    input:
        reads=lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample)
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/fastqc/fastqc.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/fastqc",
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/preassembly_fastqc_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.reads} > {log} 2>&1
        touch {output.done}
        """


rule nanoplot_long_reads:
    """
    Runs NanoPlot on ONT BAMs for long-read samples.
    """
    input:
        bam=lambda wildcards: config["input_reads"]["long_bam"].format(sample=wildcards.sample)
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/nanoplot/nanoplot.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/nanoplot",
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/preassembly_nanoplot_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} NanoPlot \
            --bam {input.bam} \
            --threads {threads} \
            --outdir {params.outdir} > {log} 2>&1
        touch {output.done}
        """


rule filtlong_long_reads:
    """
    Filters ONT reads by quality/length with Filtlong and stores a retention log.
    """
    input:
        reads=lambda wildcards: config.get("ont_reads", config["input_reads"].get("long_fastq", "")).format(sample=wildcards.sample)
    output:
        filtered=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/filtlong/filtered.fastq.gz",
        stats=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/filtlong/filtlong.stats.txt"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/filtlong",
        container=QC_CONTAINER,
        min_length=lambda wildcards: config.get("filtlong", {}).get("min_length", 1000),
        min_mean_q=lambda wildcards: config.get("filtlong", {}).get("min_mean_q", 10)
    threads: config["threads"]
    log: "logs/preassembly_filtlong_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} bash -lc "set -euo pipefail; \
            filtlong --min_length {params.min_length} --min_mean_q {params.min_mean_q} {input.reads} \
            2> {output.stats} | gzip -c > {output.filtered}" > {log} 2>&1
        """
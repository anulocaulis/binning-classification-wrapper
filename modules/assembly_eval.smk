# modules/assembly_eval.smk
# Assembly performance evaluation: metaQUAST + mapping statistics

def get_assembly_for_type(wildcards):
    assembly_inputs = config.get("assembly_inputs", {})
    if wildcards.assembly_type in assembly_inputs:
        return assembly_inputs[wildcards.assembly_type].format(sample=wildcards.sample)
    return config["input_reads"]["assembly"].format(sample=wildcards.sample)


rule metaquast_assembly_eval:
    """
    Evaluates assemblies with metaQUAST.
    """
    input:
        assembly=get_assembly_for_type
    output:
        report="{output_dir}/{sample}/assembly_eval/{assembly_type}/metaquast/report.tsv"
    params:
        outdir="{output_dir}/{sample}/assembly_eval/{assembly_type}/metaquast",
        container=QC_CONTAINER,
        extra=lambda wildcards: config.get("metaquast", {}).get("extra_args", "")
    threads: config["threads"]
    log: "logs/assembly_eval_metaquast_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} metaquast.py \
            -o {params.outdir} \
            -t {threads} \
            {params.extra} \
            {input.assembly} > {log} 2>&1
        """


rule map_short_reads_to_assembly:
    """
    Maps short reads to assemblies for mapping-based assembly performance stats.
    """
    input:
        reads=lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample),
        assembly=get_assembly_for_type
    output:
        bam="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/short_reads.sorted.bam",
        bai="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/short_reads.sorted.bam.bai"
    params:
        outdir="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping",
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/assembly_eval_map_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} bash -lc "set -euo pipefail; \
            minimap2 -t {threads} -ax sr {input.assembly} {input.reads} | \
            samtools sort -@ {threads} -o {output.bam}; \
            samtools index -@ {threads} {output.bam}" > {log} 2>&1
        """


rule samtools_mapping_stats:
    """
    Produces samtools flagstat and idxstats from mapped BAMs.
    """
    input:
        bam="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/short_reads.sorted.bam",
        bai="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/short_reads.sorted.bam.bai"
    output:
        flagstat="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/flagstat.txt",
        idxstats="{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/idxstats.txt"
    params:
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/assembly_eval_samtools_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p logs
        singularity exec {params.container} bash -lc "set -euo pipefail; \
            samtools flagstat -@ {threads} {input.bam} > {output.flagstat}; \
            samtools idxstats {input.bam} > {output.idxstats}" > {log} 2>&1
        """
# modules/assembly_eval.smk
# Assembly performance evaluation: offline QUAST + mapping statistics

def get_assembly_for_type(wildcards):
    assembly_inputs = config.get("assembly_inputs", {})
    if wildcards.assembly_type in assembly_inputs:
        return assembly_inputs[wildcards.assembly_type].format(sample=wildcards.sample)
    return config["input_reads"]["assembly"].format(sample=wildcards.sample)


def get_best_mag_reference_inputs(_wildcards):
    return config.get("best_mag_reference_inputs", [])


rule metaquast_assembly_eval:
    """
    Evaluates assemblies in offline mode with QUAST.
    """
    input:
        assembly=get_assembly_for_type
    output:
        report=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast/report.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast",
        container=QC_CONTAINER,
        extra=lambda wildcards: config.get("assembly_eval", {}).get("extra_args", "")
    threads: config["threads"]
    log: "logs/assembly_eval_metaquast_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} quast.py \
            -o {params.outdir} \
            -t {threads} \
            {params.extra} \
            {input.assembly} > {log} 2>&1
        """


rule aggregate_best_mag_references:
    """
    Aggregates curated best MAG FASTAs into a single reference file for future MetaQUAST runs.
    Not included in rule all.
    """
    input:
        mags=get_best_mag_reference_inputs
    output:
        reference=f"{OUTPUT_DIR}/references/best_mags_reference.fasta"
    log: "logs/aggregate_best_mag_references.log"
    shell:
        """
        mkdir -p $(dirname {output.reference}) logs
        if [ -z "{input.mags}" ]; then
            echo "No best_mag_reference_inputs configured in config.yaml" > {log}
            exit 1
        fi
        cat {input.mags} > {output.reference} 2> {log}
        """


rule metaquast_best_mag_reference_eval:
    """
    Evaluates assemblies with MetaQUAST against curated best MAG references.
    Not included in rule all.
    """
    input:
        assembly=get_assembly_for_type,
        reference=rules.aggregate_best_mag_references.output.reference
    output:
        report=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast_best_mags/report.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast_best_mags",
        container=QC_CONTAINER,
        extra=lambda wildcards: config.get("assembly_eval", {}).get("extra_args", "")
    threads: config["threads"]
    log: "logs/assembly_eval_metaquast_best_mags_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} metaquast.py \
            -r {input.reference} \
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
        bam=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam.bai"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping",
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
        bam=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam.bai"
    output:
        flagstat=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/flagstat.txt",
        idxstats=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/idxstats.txt"
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
# modules/nonpareil.smk
# Assembly redundancy/coverage estimation with Nonpareil.


def get_assembly_for_type(wildcards):
    assembly_inputs = config.get("assembly_inputs", {})
    if wildcards.assembly_type in assembly_inputs:
        return assembly_inputs[wildcards.assembly_type].format(sample=wildcards.sample)
    return config["input_reads"]["assembly"].format(sample=wildcards.sample)


def nonpareil_input(wildcards):
    mode = config.get("nonpareil", {}).get("input_mode", "reads")
    if mode == "assembly":
        return get_assembly_for_type(wildcards)
    return config["input_reads"]["short_interleaved"].format(sample=wildcards.sample)


rule nonpareil_assembly_eval:
    """
    Runs Nonpareil on each configured assembly type for comparison samples.
    """
    input:
        source=nonpareil_input
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/nonpareil/nonpareil.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/nonpareil",
        outprefix=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/nonpareil/nonpareil",
        staged_reads=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/nonpareil/{{sample}}.nonpareil.fastq",
        container=NONPAREIL_CONTAINER,
        algorithm=lambda wildcards: config.get("nonpareil", {}).get("algorithm", "kmer"),
        fmt=lambda wildcards: config.get("nonpareil", {}).get("format", "fastq"),
        extra=lambda wildcards: config.get("nonpareil", {}).get("extra_args", "")
    threads: 1
    log: "logs/assembly_eval_nonpareil_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        NP_INPUT={input.source}
        if [ "${{NP_INPUT##*.}}" = "gz" ]; then
            if [ ! -s {params.staged_reads} ]; then
                gzip -dc "$NP_INPUT" > {params.staged_reads}
            fi
            NP_INPUT={params.staged_reads}
        fi
        singularity exec {params.container} sh -c "set -euo pipefail; \
            command -v nonpareil >/dev/null; \
            nonpareil \
                -s $NP_INPUT \
                -T {params.algorithm} \
                -f {params.fmt} \
                -b {params.outprefix} \
                {params.extra}" > {log} 2>&1
        touch {output.done}
        """

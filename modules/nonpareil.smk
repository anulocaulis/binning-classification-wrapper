# modules/nonpareil.smk
# Assembly redundancy/coverage estimation with Nonpareil.


rule nonpareil_assembly_eval:
    """
    Runs Nonpareil for assembly redundancy/coverage estimation.
    
    **Inputs:**
    - source: Trimmed and filtered interleaved FASTQ from readmap_prep module.
              Contains polyG-filtered, quality-trimmed reads for accurate coverage analysis.
    
    **Outputs:**
    - done: Marker file indicating successful Nonpareil completion.
            The tool also produces *.npo, *.npe, and *.npc files in the outdir.
    
    **Logic:**
    Nonpareil estimates how much of the sample diversity has been sequenced based on
    read redundancy patterns. If input is gzip-compressed, it's decompressed to a temporary
    FASTQ file for Nonpareil processing. The algorithm uses k-mer overlap analysis
    (kmer mode by default) to determine coverage/redundancy without requiring alignment.
    """
    input:
        source=f"trimmed_reads/{{sample}}_interleaved_trimmed_polyG_filtered.fastq.gz"
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

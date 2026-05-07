# ---
# Step 1: Trim Illumina reads after fastqc
# ---
rule trim_illumina:
    """
    Trims Illumina reads for adapters and quality using bbduk.sh from BBTools.
    Takes interleaved FASTQ as input and outputs a trimmed interleaved FASTQ.
    """
    input:
        reads=config["input_reads"]["short_interleaved"]
    output:
        trimmed_reads="trimmed_reads/{sample}_interleaved_trimmed.fastq.gz"
    container:
        QC_CONTAINER
    params:
        container_path = QC_CONTAINER,
        # - ref=adapters: Uses built-in Illumina adapter list
        # - ktrim=r k=23 mink=11 hdist=1: Standard right-end adapter trimming
        # - tpe: Trim Paired Ends
        # - tbo: Trim By Overlap
        # - qtrim=rl trimq=10: Quality trimming from left and right ends
        # - minlen=50: Discards pairs where either read is < 50 bp after trimming
        bbduk_opts="ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 minlen=50"
    threads:
        config["threads"]
    log:
        "logs/trim_illumina_{sample}.log"
    shell:
        """
        WORKDIR=$(pwd)
        singularity exec -B "$WORKDIR:$WORKDIR" -W "$WORKDIR" {params.container_path} bbduk.sh in={input.reads} out={output.trimmed_reads} threads={threads} {params.bbduk_opts} 2>> {log}
        """

# ---
# Step 2: Remove polyG tails from Illumina reads
# ---
rule trim_polyg:
    """
    Removes polyG tails from Illumina reads using BBTools polyfilter.sh
    """
    input:
        reads = "trimmed_reads/{sample}_interleaved_trimmed.fastq.gz"
    output:
        filtered_reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    container:
        QC_CONTAINER
    params:
        container_path = QC_CONTAINER
    threads:
        config["threads"]
    log:
        "logs/trim_polyg_{sample}.log"
    shell:
        """
        WORKDIR=$(pwd)
        singularity exec -B "$WORKDIR:$WORKDIR" {params.container_path} polyfilter.sh in={input.reads} out={output.filtered_reads} 2>> {log}
        """

# ---
# Step 3: de-interleave reads after trimming
# ---
rule deinterleave_illumina:
    """De-interleaves trimmed Illumina reads into separate R1 and R2 FASTQ files using BBTools reformat.sh"""
    input:
        reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        r1 = "trimmed_reads/{sample}/{sample}_R1.fastq.gz",
        r2 = "trimmed_reads/{sample}/{sample}_R2.fastq.gz"
    container:
        QC_CONTAINER
    params:
        container_path = QC_CONTAINER
    threads:
        config["threads"]
    log:
        "logs/deinterleave_{sample}.log"
    shell:
        """
        WORKDIR=$(pwd)
        mkdir -p trimmed_reads/{wildcards.sample} logs
        singularity exec -B "$WORKDIR:$WORKDIR" -W "$WORKDIR" {params.container_path} reformat.sh in={input.reads} out1={output.r1} out2={output.r2} threads={threads} 2>> {log}
        """


def resolve_assembly_path(sample, assembly_type):
    """
    Resolve the actual assembly file path for a given sample and assembly type.
    Tries multiple naming conventions since different assemblers use different names.
    """
    import os

    assembly_root = config.get("assembly_root_dir", "/storage/biology/projects/miller-lowry/beitner/data/assemblies")
    assembly_short = assembly_type.replace("assembly.", "")
    sample_lower = sample.lower()

    canonical_base = f"{assembly_root}/{assembly_short}"
    canonical_candidates = [
        f"{canonical_base}/{sample}.assembly.fasta",
        f"{canonical_base}/{sample_lower}.assembly.fasta",
        f"{canonical_base}/{sample}.assembly.fa",
        f"{canonical_base}/{sample_lower}.assembly.fa",
        f"{canonical_base}/{sample}.{assembly_short}.fasta",
        f"{canonical_base}/{sample_lower}.{assembly_short}.fasta",
        f"{canonical_base}/{sample}.{assembly_short}.fa",
        f"{canonical_base}/{sample_lower}.{assembly_short}.fa",
    ]

    legacy_base = f"{assembly_root}/{sample}/{assembly_type}"
    legacy_candidates = [
        f"{legacy_base}/assembly.fasta",
        f"{legacy_base}/contigs.fasta",
        f"{legacy_base}/final.contigs.fa",
        f"{legacy_base}/{sample}_polished.fasta",
        f"{legacy_base}/{assembly_short}.contigs.fasta",
        f"{legacy_base}/contigs.ge1000.fa",
    ]

    candidates = canonical_candidates + legacy_candidates
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(f"No assembly found for {sample}/{assembly_type}. Tried: {candidates}")


def configured_cobinning_samples(assembly_subdir):
    """Return configured sample list for an assembly_subdir, preserving order."""
    assembly_map = config.get("cobinning_samples_by_assembly", {})
    requested = assembly_map.get(assembly_subdir)
    if requested:
        return [s for s in requested if s in ALL_SAMPLES]
    return ALL_SAMPLES


def samples_with_assembly(assembly_subdir):
    """Return samples that have a resolvable assembly for the given assembly_subdir."""
    resolved = []
    for sample in configured_cobinning_samples(assembly_subdir):
        try:
            resolve_assembly_path(sample, assembly_subdir)
            resolved.append(sample)
        except FileNotFoundError:
            continue
    if not resolved:
        raise ValueError(f"No samples have a resolvable assembly for {assembly_subdir}")
    return resolved


def cobinning_staged_read_pairs(wildcards):
    """
    Return the staged read pairs for all samples in the cobinning group for the given assembly_subdir.
    This is used as input for the MetaWRAP co-binning step, which requires all reads to be specified together.
    """
    return " ".join(
        f"{OUTPUT_DIR}/{sample}/binning/work/{sample}_1.fastq::{OUTPUT_DIR}/{sample}/binning/work/{sample}_2.fastq"
        for sample in samples_with_assembly(wildcards.assembly_subdir)
    )


def cobinning_runtime_minutes(assembly_subdir):
    per_assembly = config.get("cobinning_runtime_min_by_assembly", {})
    if assembly_subdir in per_assembly:
        return int(per_assembly[assembly_subdir])
    return int(config.get("cobinning_runtime_min", 2880))


def cobinning_threads(assembly_subdir):
    per_assembly = config.get("cobinning_threads_by_assembly", {})
    if assembly_subdir in per_assembly:
        return int(per_assembly[assembly_subdir])
    return int(config.get("cobinning_threads", config["threads"]))


def cobinning_mem_mb(assembly_subdir):
    per_assembly = config.get("cobinning_mem_mb_by_assembly", {})
    if assembly_subdir in per_assembly:
        return int(per_assembly[assembly_subdir])
    return int(config.get("cobinning_mem_mb", 360000))


# ---
# Shared prep: stage split reads for MetaWRAP naming requirements
# ---
rule stage_reads_for_metawrap:
    """
    MetaWRAP binning requires uncompressed reads named *_1.fastq and *_2.fastq.
    This rule decompresses the preassembly_qc split reads and stages them.
    """
    input:
        r1=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/split_reads/{{sample}}_R1.fastq.gz",
        r2=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/split_reads/{{sample}}_R2.fastq.gz"
    output:
        r1=f"{OUTPUT_DIR}/{{sample}}/binning/work/{{sample}}_1.fastq",
        r2=f"{OUTPUT_DIR}/{{sample}}/binning/work/{{sample}}_2.fastq"
    log: "logs/stage_reads_metawrap_{sample}.log"
    shell:
        """
        : > {log}
        mkdir -p $(dirname {output.r1}) logs
        gzip -dc {input.r1} > {output.r1} 2>> {log}
        gzip -dc {input.r2} > {output.r2} 2>> {log}
        """


# modules/comingled_binning.smk
# Multi-split (sample-wise assembly, multi-sample coverage) binning.
# Concatenates per-sample assemblies for a given assembly type, prefixes contigs by sample,
# and runs pooled MetaBAT2/CONCOCT (via MetaWRAP) plus pooled VAMB.


def comingled_assembly_path(sample, assembly_type):
    assembly_inputs = config.get("assembly_inputs", {})
    if assembly_type in assembly_inputs:
        return assembly_inputs[assembly_type].format(sample=sample)
    return config["input_reads"]["assembly"].format(sample=sample)


def comingled_assemblies(wildcards):
    return [comingled_assembly_path(sample, wildcards.assembly_type) for sample in ALL_SAMPLES]


def comingled_sample_assembly_pairs(wildcards):
    return " ".join(
        f"{sample}::{comingled_assembly_path(sample, wildcards.assembly_type)}" for sample in ALL_SAMPLES
    )


rule comingled_concat_assemblies:
    """
    Concatenate all sample assemblies for one assembly type into a shared FASTA.
    Each contig header is prefixed with sample ID to avoid ID collisions.
    """
    input:
        assemblies=comingled_assemblies
    output:
        fasta=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/coassembly/combined.prefixed.fa"
    params:
        sample_assemblies=comingled_sample_assembly_pairs
    log: "logs/comingled_concat_{assembly_type}.log"
    shell:
        """
        mkdir -p $(dirname {output.fasta}) logs
        : > {output.fasta}
        : > {log}

        for pair in {params.sample_assemblies}; do
            sample="${{pair%%::*}}"
            fasta="${{pair#*::}}"
            if [ ! -s "$fasta" ]; then
                echo "Missing/empty assembly for $sample: $fasta" >> {log}
                exit 1
            fi
            awk -v pref="$sample|" '
                /^>/ {
                    h = substr($0, 2)
                    sub(/[ \t].*$/, "", h)
                    print ">" pref h
                    next
                }
                { print }
            ' "$fasta" >> {output.fasta}
        done

        if [ ! -s {output.fasta} ]; then
            echo "Combined assembly is empty: {output.fasta}" >> {log}
            exit 1
        fi
        """


rule comingled_metawrap_binning:
    """
    Run pooled MetaWRAP binning (MetaBAT2 + CONCOCT) against the shared assembly
    using interleaved reads from all samples.
    """
    input:
        fasta=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/coassembly/combined.prefixed.fa",
        reads=expand(config["input_reads"]["short_interleaved"], sample=ALL_SAMPLES)
    output:
        metabat2=directory(f"{OUTPUT_DIR}/comingled/{{assembly_type}}/metawrap/metabat2_bins"),
        concoct=directory(f"{OUTPUT_DIR}/comingled/{{assembly_type}}/metawrap/concoct_bins"),
        done=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/metawrap/comingled_metawrap.done"
    params:
        outdir=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/metawrap",
        staged_reads=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/metawrap/.staged_reads",
        container=METAWRAP_CONTAINER,
        metawrap_exec=lambda wildcards: config.get("metawrap_executable", "/usr/local/bin/metawrap"),
        min_contig_len=lambda wildcards: str(config.get("comingled_binning", {}).get("min_contig_len", 1000))
    threads: config["threads"]
    log: "logs/comingled_metawrap_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} {params.staged_reads} logs
        : > {log}

        READS_ARGS=""
        for r in {input.reads}; do
            base=$(basename "$r")
            staged="{params.staged_reads}/${base%.gz}"
            if [ "${r##*.}" = "gz" ]; then
                if [ ! -s "$staged" ]; then
                    gzip -dc "$r" > "$staged"
                fi
                READS_ARGS="$READS_ARGS $staged"
            else
                READS_ARGS="$READS_ARGS $r"
            fi
        done

        singularity exec {params.container} {params.metawrap_exec} binning \
            -o {params.outdir} \
            -t {threads} \
            -a {input.fasta} \
            --metabat2 --concoct \
            --interleaved \
            -l {params.min_contig_len} \
            $READS_ARGS >> {log} 2>&1

        if [ ! -d {output.metabat2} ] || [ ! -d {output.concoct} ]; then
            echo "Expected pooled metabat2/concoct output dirs were not created." >> {log}
            exit 1
        fi

        touch {output.done}
        """


rule comingled_map_sample_to_combined:
    """
    Map each sample's reads to the shared combined assembly for pooled VAMB abundance.
    """
    input:
        fasta=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/coassembly/combined.prefixed.fa",
        reads=lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample)
    output:
        bam=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/mapping/{{sample}}.sorted.bam.bai"
    params:
        outdir=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/mapping",
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/comingled_map_{assembly_type}_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} sh -c "set -euo pipefail; \
            export PATH=/opt/conda/envs/qc_env/bin:/opt/conda/envs/base_tools/bin:/opt/conda/bin:$PATH; \
            minimap2 -t {threads} -ax sr {input.fasta} {input.reads} | \
            samtools sort -@ {threads} -o {output.bam}; \
            samtools index -@ {threads} {output.bam}" > {log} 2>&1
        """


rule comingled_vamb_abundance_tsv:
    """
    Build a combined VAMB abundance TSV with one abundance column per sample.
    """
    input:
        fasta=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/coassembly/combined.prefixed.fa",
        bams=lambda wildcards: expand(
            f"{OUTPUT_DIR}/comingled/{{assembly_type}}/mapping/{{sample}}.sorted.bam",
            assembly_type=wildcards.assembly_type,
            sample=ALL_SAMPLES,
        ),
        bais=lambda wildcards: expand(
            f"{OUTPUT_DIR}/comingled/{{assembly_type}}/mapping/{{sample}}.sorted.bam.bai",
            assembly_type=wildcards.assembly_type,
            sample=ALL_SAMPLES,
        )
    output:
        tsv=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb/combined.abundance.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb",
        container=QC_CONTAINER,
        samples=" ".join(ALL_SAMPLES)
    threads: config["threads"]
    log: "logs/comingled_vamb_abundance_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        tmpdir=$(mktemp -d)

        awk '/^>/ {{h=substr($0,2); sub(/[ \t].*$/, "", h); print h}}' {input.fasta} | sort -u > "$tmpdir/contigs.tsv"
        cp "$tmpdir/contigs.tsv" "$tmpdir/table.tsv"

        for s in {params.samples}; do
            bam="{OUTPUT_DIR}/comingled/{wildcards.assembly_type}/mapping/${{s}}.sorted.bam"
            singularity exec {params.container} samtools idxstats "$bam" | \
                awk 'BEGIN{{OFS="\t"}} $1 != "*" {{len=$2+0; mapped=$3+0; abund=(len>0 ? mapped/len : 0); print $1, abund}}' | \
                sort -k1,1 > "$tmpdir/${{s}}.tsv"
            join -t $'\t' -a 1 -e 0 -o 1.1,2.2 "$tmpdir/table.tsv" "$tmpdir/${{s}}.tsv" > "$tmpdir/table.next.tsv"
            mv "$tmpdir/table.next.tsv" "$tmpdir/table.tsv"
        done

        printf 'contigname\t{params.samples}\n' | sed 's/ /\t/g' > {output.tsv}
        cat "$tmpdir/table.tsv" >> {output.tsv}

        if [ $(wc -l < {output.tsv}) -le 1 ]; then
            echo "Combined abundance TSV has no contig rows: {output.tsv}" >&2
            rm -rf "$tmpdir"
            exit 1
        fi

        rm -rf "$tmpdir"
        """


rule comingled_vamb_bin_default:
    """
    Run pooled VAMB binning using combined assembly and multi-sample abundance matrix.
    """
    input:
        fasta=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/coassembly/combined.prefixed.fa",
        abundance=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb/combined.abundance.tsv"
    output:
        bins=directory(f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb/vamb_bins"),
        done=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb/comingled_vamb.done"
    params:
        outdir=f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb/vamb_bins",
        container=VAMB_CONTAINER,
        min_contig_len=lambda wildcards: str(config.get("comingled_binning", {}).get("vamb_min_contig_len", 2000)),
        threads=lambda wildcards: str(config.get("comingled_binning", {}).get("vamb_threads", config.get("threads", 8))),
        norefcheck=lambda wildcards: bool(config.get("comingled_binning", {}).get("vamb_norefcheck", False))
    threads: config["threads"]
    log: "logs/comingled_vamb_bin_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        norefcheck_flag=""
        if [ "{params.norefcheck}" = "True" ]; then
            norefcheck_flag="--norefcheck"
        fi

        singularity exec {params.container} vamb bin default \
            --outdir {params.outdir} \
            --fasta {input.fasta} \
            --abundance_tsv {input.abundance} \
            -m {params.min_contig_len} \
            -p {params.threads} \
            $norefcheck_flag > {log} 2>&1

        if [ ! -s {params.outdir}/clusters.tsv ]; then
            echo "Pooled VAMB finished but clusters.tsv is missing in {params.outdir}" >&2
            exit 1
        fi

        touch {output.done}
        """

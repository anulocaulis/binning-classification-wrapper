# modules/multiqc.smk
# Aggregate all pipeline outputs into a consolidated MultiQC report.

MULTIQC_INPUTS_DIR = f"{OUTPUT_DIR}/multiqc_inputs"


rule multiqc_collect:
    """
    Symlinks all text/TSV/CSV analysis outputs into a per-sample staging directory
    (multiqc_inputs/) so MultiQC scans a clean, purpose-built tree rather than
    the entire binning_outputs/ hierarchy.
    """
    input:
        BASE_TARGETS
    output:
        done=f"{MULTIQC_INPUTS_DIR}/collected.done"
    params:
        outdir=MULTIQC_INPUTS_DIR,
        assembly_types=" ".join(ASSEMBLY_TYPES),
        compare_samples=" ".join(COMPARE_SAMPLES),
        all_samples=" ".join(ALL_SAMPLES),
    log: "logs/multiqc_collect.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        for sample in {params.all_samples}; do
            sdir="{params.outdir}/$sample"
            mkdir -p "$sdir"

            # FastQC (directory symlink)
            src="{OUTPUT_DIR}/$sample/preassembly_qc/fastqc"
            [ -d "$src" ] && ln -sfn "$(realpath "$src")" "$sdir/fastqc" || true

            # NanoPlot (long-read samples only)
            src="{OUTPUT_DIR}/$sample/preassembly_qc/nanoplot"
            [ -d "$src" ] && ln -sfn "$(realpath "$src")" "$sdir/nanoplot" || true

            # Filtlong stats
            src="{OUTPUT_DIR}/$sample/preassembly_qc/filtlong/filtlong.stats.txt"
            [ -f "$src" ] && ln -sf "$(realpath "$src")" "$sdir/filtlong_stats.txt" || true

            # Kraken2 report
            src="{OUTPUT_DIR}/$sample/classification/kraken2/report.txt"
            [ -f "$src" ] && ln -sf "$(realpath "$src")" "$sdir/kraken2_report.txt" || true

            # GTDB-Tk summary
            src="{OUTPUT_DIR}/$sample/classification/gtdbtk/gtdbtk.bac120.summary.tsv"
            [ -f "$src" ] && ln -sf "$(realpath "$src")" "$sdir/gtdbtk_summary.tsv" || true

            # CheckM2 quality
            src="{OUTPUT_DIR}/$sample/checkm2/quality_report.tsv"
            [ -f "$src" ] && ln -sf "$(realpath "$src")" "$sdir/checkm2_quality_report.tsv" || true

            # BUSCO outputs (directory symlink)
            src="{OUTPUT_DIR}/$sample/mag_integrity/busco"
            [ -d "$src" ] && ln -sfn "$(realpath "$src")" "$sdir/busco" || true

            # GUNC outputs (directory symlink)
            src="{OUTPUT_DIR}/$sample/mag_integrity/gunc"
            [ -d "$src" ] && ln -sfn "$(realpath "$src")" "$sdir/gunc" || true

            # Prokka outputs (directory symlink)
            src="{OUTPUT_DIR}/$sample/functional/prokka"
            [ -d "$src" ] && ln -sfn "$(realpath "$src")" "$sdir/prokka" || true
        done

        # Per-assembly-type mapping and metaquast outputs
        for sample in {params.compare_samples}; do
            sdir="{params.outdir}/$sample"
            mkdir -p "$sdir"
            for atype in {params.assembly_types}; do
                adir="$sdir/assembly_eval_$atype"
                mkdir -p "$adir"

                src="{OUTPUT_DIR}/$sample/assembly_eval/$atype/metaquast/report.tsv"
                [ -f "$src" ] && ln -sf "$(realpath "$src")" "$adir/metaquast_report.tsv" || true

                src="{OUTPUT_DIR}/$sample/assembly_eval/$atype/mapping/flagstat.txt"
                [ -f "$src" ] && ln -sf "$(realpath "$src")" "$adir/flagstat.txt" || true

                src="{OUTPUT_DIR}/$sample/assembly_eval/$atype/mapping/idxstats.txt"
                [ -f "$src" ] && ln -sf "$(realpath "$src")" "$adir/idxstats.txt" || true
            done
        done

        echo "multiqc_collect done: $(date)" >> {log}
        touch {output.done}
        """


rule multiqc_report:
    """
    Runs MultiQC across the staged multiqc_inputs/ directory.
    """
    input:
        collected=f"{MULTIQC_INPUTS_DIR}/collected.done"
    output:
        report=f"{OUTPUT_DIR}/multiqc/multiqc_report.html",
        data_dir=directory(f"{OUTPUT_DIR}/multiqc/multiqc_data")
    params:
        outdir=f"{OUTPUT_DIR}/multiqc",
        inputs_dir=MULTIQC_INPUTS_DIR,
        container=QC_CONTAINER,
        title=config.get("multiqc", {}).get("title", "Binning Classification Workflow"),
        comment=config.get("multiqc", {}).get("comment", "Consolidated QC and analysis summary"),
        extra_args=config.get("multiqc", {}).get("extra_args", "")
    log:
        "logs/multiqc.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} sh -c "set -euo pipefail; \
            export PATH=/opt/conda/bin:$PATH; \
            command -v multiqc >/dev/null; \
            multiqc {params.inputs_dir} logs \
                -o {params.outdir} \
                -n multiqc_report.html \
                --force \
                --title '{params.title}' \
                --comment '{params.comment}' \
                {params.extra_args}" > {log} 2>&1
        """

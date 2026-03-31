# modules/multiqc.smk
# Aggregate all pipeline outputs into a consolidated MultiQC report.

rule multiqc_report:
    """
    Runs MultiQC across workflow outputs/logs after all primary targets complete.
    """
    input:
        BASE_TARGETS
    output:
        report=f"{OUTPUT_DIR}/multiqc/multiqc_report.html",
        data_dir=directory(f"{OUTPUT_DIR}/multiqc/multiqc_data")
    params:
        outdir=f"{OUTPUT_DIR}/multiqc",
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
            multiqc {OUTPUT_DIR} logs \
                -o {params.outdir} \
                -n multiqc_report.html \
                --force \
                --title '{params.title}' \
                --comment '{params.comment}' \
                {params.extra_args}" > {log} 2>&1
        """

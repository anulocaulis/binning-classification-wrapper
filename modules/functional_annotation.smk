# modules/functional_annotation.smk
# Functional annotation of refined bins

rule prokka_annotate_bins:
    """
    Annotates refined bins with Prokka to estimate functional gene content.
    """
    input:
        bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/functional/prokka/prokka.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/functional/prokka",
        container=CLASSIFICATION_CONTAINER,
        kingdom=lambda wildcards: config.get("prokka", {}).get("kingdom", "Bacteria")
    threads: config["threads"]
    log: "logs/functional_prokka_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} bash -lc "set -euo pipefail; \
            for f in {input.bins_dir}/*.fa; do \
                name=$(basename \"$f\" .fa); \
                prokka --outdir {params.outdir}/$name --prefix $name \
                    --cpus {threads} --kingdom {params.kingdom} \"$f\"; \
            done" > {log} 2>&1
        touch {output.done}
        """
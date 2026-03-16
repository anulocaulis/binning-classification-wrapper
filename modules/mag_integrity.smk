# modules/mag_integrity.smk
# Biological integrity checks for refined MAG bins

rule gunc_assess_bins:
    """
    Detects potential chimerism/contamination in refined bins with GUNC.
    """
    input:
        bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/gunc/gunc.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/gunc",
        container=CLASSIFICATION_CONTAINER,
        db_arg=lambda wildcards: ("--db_file " + config.get("gunc_db", "")) if config.get("gunc_db", "") else ""
    threads: config["threads"]
    log: "logs/mag_integrity_gunc_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} bash -lc "set -euo pipefail; \
            gunc run --input_dir {input.bins_dir} --out_dir {params.outdir} --threads {threads} \
            {params.db_arg}" > {log} 2>&1
        touch {output.done}
        """


rule busco_assess_bins:
    """
    Runs BUSCO on refined bins as a secondary conserved-gene completeness check.
    """
    input:
        bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/busco/busco.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/busco",
        container=CLASSIFICATION_CONTAINER,
        lineage=lambda wildcards: config.get("busco", {}).get("lineage", "bacteria_odb10"),
        mode=lambda wildcards: config.get("busco", {}).get("mode", "genome")
    threads: config["threads"]
    log: "logs/mag_integrity_busco_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} bash -lc "set -euo pipefail; \
            for f in {input.bins_dir}/*.fa; do \
                name=$(basename \"$f\" .fa); \
                busco -i \"$f\" -o \"$name\" -m {params.mode} -l {params.lineage} \
                    --cpu {threads} --out_path {params.outdir}; \
            done" > {log} 2>&1
        touch {output.done}
        """
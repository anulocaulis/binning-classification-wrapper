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
        container=CLASSIFICATION_GTD_GUNC_CONTAINER,
        executable="gunc",
        db_arg=lambda wildcards: ("--db_file " + config.get("gunc_db", "")) if config.get("gunc_db", "") else ""
    threads: config["threads"]
    log: "logs/mag_integrity_gunc_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        if [ -z "{params.db_arg}" ]; then
            echo "Skipping GUNC because gunc_db is not configured." > {log}
            touch {output.done}
            exit 0
        fi
        singularity exec {params.container} sh -c "command -v {params.executable} >/dev/null 2>&1" > {log} 2>&1
        singularity exec {params.container} {params.executable} run \
            --input_dir {input.bins_dir} --out_dir {params.outdir} --threads {threads} \
            {params.db_arg} >> {log} 2>&1
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
        executable="/opt/conda/envs/base_tools/bin/busco",
        lineage=lambda wildcards: config.get("busco", {}).get("lineage", "bacteria_odb10"),
        mode=lambda wildcards: config.get("busco", {}).get("mode", "genome")
    threads: config["threads"]
    log: "logs/mag_integrity_busco_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}
        for f in {input.bins_dir}/*.fa; do
            name=$(basename "$f" .fa)
            singularity exec {params.container} {params.executable} \
                -i "$f" -o "$name" -m {params.mode} -l {params.lineage} \
                --cpu {threads} --out_path {params.outdir} >> {log} 2>&1
        done
        touch {output.done}
        """
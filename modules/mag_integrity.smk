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
            summary="{params.outdir}/$name/short_summary.specific.{params.lineage}.$name.txt"
            if [ -s "$summary" ]; then
                echo "Skipping completed BUSCO run for $name" >> {log}
                continue
            fi
            rm -rf "{params.outdir}/$name"
            singularity exec {params.container} {params.executable} \
                -i "$f" -o "$name" -m {params.mode} -l {params.lineage} -f \
                --cpu {threads} --out_path {params.outdir} >> {log} 2>&1
        done
        touch {output.done}
        """


rule magqual_label_mimag:
    """
    Runs Bakta on refined bins inside magqual.sif.
    This intentionally skips MAGqual's internal CheckM1 workflow.
    """
    input:
        bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins",
        assembly=lambda wildcards: config["input_reads"]["assembly"].format(sample=wildcards.sample)
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/magqual/mimag.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/magqual",
        container=MAGQUAL_CONTAINER,
        bakta_db=lambda wildcards: config.get("magqual", {}).get("bakta_db", "")
    threads: config["threads"]
    log: "logs/mag_integrity_magqual_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        if [ -z "{params.bakta_db}" ]; then
            echo "config.magqual.bakta_db must be set for Bakta-only mode." >> {log}
            exit 1
        fi
        if [ ! -d "{params.bakta_db}" ]; then
            echo "Configured config.magqual.bakta_db path does not exist: {params.bakta_db}" >> {log}
            exit 1
        fi

        container_abs=$(realpath {params.container})
        bins_abs=$(realpath {input.bins_dir})
        bakta_db_abs=$(realpath {params.bakta_db})
        out_abs=$(realpath {params.outdir})

        mkdir -p "$out_abs/analysis/bakta/{wildcards.sample}"

        # bakta 1.8.1 in magqual.sif uses pyrodigal.OrfFinder (pre-3.x API).
        # Install pyrodigal==2.1.0 to a local dir and override via PYTHONPATH.
        pyrodigal_fix="$out_abs/.pyrodigal_fix"
        mkdir -p "$pyrodigal_fix"
        singularity exec \
            -B "$out_abs:$out_abs" \
            "$container_abs" \
            sh -lc "PATH=/opt/conda/envs/magqual/bin:$PATH; pip install 'pyrodigal==2.1.0' -q --target '$pyrodigal_fix'" >> {log} 2>&1

        for f in "$bins_abs"/*.fa; do
            name=$(basename "$f" .fa)
            bin_out="$out_abs/analysis/bakta/{wildcards.sample}/$name"

            if [ -s "$bin_out/$name.tsv" ] && [ -s "$bin_out/$name.txt" ] && [ -s "$bin_out/$name.faa" ]; then
                echo "Skipping completed Bakta run for $name" >> {log}
                continue
            fi

            rm -rf "$bin_out"

            singularity exec \
                -B "$bins_abs:$bins_abs" \
                -B "$out_abs:$out_abs" \
                -B "$bakta_db_abs:$bakta_db_abs" \
                "$container_abs" \
                sh -lc "PATH=/opt/conda/envs/magqual/bin:$PATH; PYTHONPATH='$pyrodigal_fix' bakta --force --db '$bakta_db_abs' --output '$bin_out' --prefix '$name' --threads {threads} '$f'" >> {log} 2>&1
        done

        touch {output.done}
        """
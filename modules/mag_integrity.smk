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


rule magqual_label_mimag:
    """
    Runs MAGqual independently to produce native MIMAG labeling for refined bins.
    """
    input:
        bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins",
        assembly=lambda wildcards: config["input_reads"]["assembly"].format(sample=wildcards.sample)
    output:
        labels=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/magqual/analysis/genome_bins/{{sample}}_qual_MAGs.txt",
        done=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/magqual/mimag.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/mag_integrity/magqual",
        container=MAGQUAL_CONTAINER,
        bakta_db=lambda wildcards: config.get("magqual", {}).get("bakta_db", ""),
        checkm_db=lambda wildcards: config.get("magqual", {}).get("checkm_db", "")
    threads: config["threads"]
    log: "logs/mag_integrity_magqual_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        container_abs=$(realpath {params.container})
        log_abs=$(realpath {log})

        assembly_abs=$(realpath {input.assembly})
        bins_abs=$(realpath {input.bins_dir})

        bakta_arg=""
        bakta_bind=""
        if [ -n "{params.bakta_db}" ]; then
            if [ ! -d "{params.bakta_db}" ]; then
                echo "Configured config.magqual.bakta_db path does not exist: {params.bakta_db}" >> {log}
                exit 1
            fi
            bakta_arg="--baktadb {params.bakta_db}"
            bakta_bind="-B {params.bakta_db}:{params.bakta_db}"
        else
            echo "config.magqual.bakta_db is empty; MAGqual will download/use its own Bakta light database in the MAGqual run directory." >> {log}
        fi

        checkm_arg=""
        if [ -n "{params.checkm_db}" ]; then
            if [ ! -d "{params.checkm_db}" ]; then
                echo "Configured config.magqual.checkm_db path does not exist: {params.checkm_db}" >> {log}
                exit 1
            fi
            checkm_arg="--checkmdb {params.checkm_db}"
        fi

        cd {params.outdir}
        singularity exec -B "$PWD:$PWD" $bakta_bind -B "$bins_abs:$bins_abs" -B "$assembly_abs:$assembly_abs" "$container_abs" \
            sh -c "magqual -a '$assembly_abs' -b '$bins_abs' -p '{wildcards.sample}' -j {threads} $bakta_arg $checkm_arg" >> "$log_abs" 2>&1

        test -s analysis/genome_bins/{wildcards.sample}_qual_MAGs.txt
        touch mimag.done
        """
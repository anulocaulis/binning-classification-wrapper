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
        executable="/opt/conda/envs/base_tools/bin/prokka",
        kingdom=lambda wildcards: config.get("prokka", {}).get("kingdom", "Bacteria")
    threads: config["threads"]
    log: "logs/functional_prokka_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}
        for f in {input.bins_dir}/*.fa; do
            name=$(basename "$f" .fa)
            outdir="{params.outdir}/$name"
            gff="$outdir/$name.gff"

            if [[ -s "$gff" ]]; then
                echo "[$(date '+%F %T')] Skipping $name (existing annotation found: $gff)" >> {log}
                continue
            fi

            rm -rf "$outdir"
            singularity exec {params.container} {params.executable} \
                --outdir "$outdir" --prefix "$name" \
                --cpus {threads} --kingdom {params.kingdom} "$f" >> {log} 2>&1
        done
        touch {output.done}
        """
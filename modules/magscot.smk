# modules/magscot.smk
# MAGScoT bin scoring/refinement module

rule magscot_generate_contig_bin_map:
    """
    Generate MAGScoT -i input table with three columns (no header):
    bin_name, contig_id, bin_set_label.
    """
    input:
        metabat2=f"{OUTPUT_DIR}/{{sample}}/binning/metabat2_bins",
        maxbin2=f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2_bins",
        concoct=f"{OUTPUT_DIR}/{{sample}}/binning/concoct_bins"
    output:
        table=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/{{sample}}.input.tsv",
        done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/contig_bin_map.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot"
    log: "logs/magscot_contig_bin_map_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}
        : > {output.table}

        for f in {input.metabat2}/*.fa {input.metabat2}/*.fasta; do
            [ -e "$f" ] || continue
            b=$(basename "$f")
            awk -v bin="$b" -v set="metabat2" '/^>/{{c=substr($0,2); sub(/[ \t].*$/, "", c); print bin "\t" c "\t" set}}' "$f" >> {output.table}
        done
        for f in {input.maxbin2}/*.fa {input.maxbin2}/*.fasta; do
            [ -e "$f" ] || continue
            b=$(basename "$f")
            awk -v bin="$b" -v set="maxbin2" '/^>/{{c=substr($0,2); sub(/[ \t].*$/, "", c); print bin "\t" c "\t" set}}' "$f" >> {output.table}
        done
        for f in {input.concoct}/*.fa {input.concoct}/*.fasta; do
            [ -e "$f" ] || continue
            b=$(basename "$f")
            awk -v bin="$b" -v set="concoct" '/^>/{{c=substr($0,2); sub(/[ \t].*$/, "", c); print bin "\t" c "\t" set}}' "$f" >> {output.table}
        done

        if [ ! -s {output.table} ]; then
            echo "Generated MAGScoT contig-to-bin table is empty: {output.table}" >> {log}
            exit 1
        fi

        touch {output.done}
        """


rule magscot_predict_proteins:
    """
    Predict proteins from assembly contigs using prodigal in metagenome mode (-p meta).
    """
    input:
        assembly=lambda wildcards: config["input_reads"]["assembly"].format(sample=wildcards.sample)
    output:
        faa=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/{{sample}}.proteins.faa",
        done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/prodigal.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot",
        container="containers/magscot.sif"
    log: "logs/magscot_prodigal_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        singularity exec {params.container} prodigal \
            -i {input.assembly} \
            -a {output.faa} \
            -p meta >> {log} 2>&1

        if [ ! -s {output.faa} ]; then
            echo "Prodigal produced no proteins: {output.faa}" >> {log}
            exit 1
        fi

        touch {output.done}
        """


rule magscot_generate_markers:
    """
    Generate MAGScoT --hmm marker mapping TSV (gene_id, marker, e-value)
    from predicted proteins using hmmsearch.
    """
    input:
        faa=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/{{sample}}.proteins.faa",
        prodigal_done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/prodigal.done"
    output:
        hmm=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/{{sample}}.markers.tsv",
        done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/markers.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot",
        container="containers/magscot.sif",
        marker_hmm_db=lambda wildcards: config.get("magscot", {}).get("marker_hmm_db", "") or config.get("magscot", {}).get("hmm", ""),
        hmmer_evalue=lambda wildcards: str(config.get("magscot", {}).get("hmmer_evalue", 1e-10))
    threads: config["threads"]
    log: "logs/magscot_markers_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        if [ -z "{params.marker_hmm_db}" ]; then
            echo "config.magscot.marker_hmm_db (or legacy config.magscot.hmm) is required for marker generation." >> {log}
            exit 1
        fi
        if [ ! -e "{params.marker_hmm_db}" ]; then
            echo "MAGScoT marker HMM DB not found: {params.marker_hmm_db}" >> {log}
            exit 1
        fi

        tmpdir="{params.outdir}/work_markers"
        rm -rf "$tmpdir"
        mkdir -p "$tmpdir"
        tbl="$tmpdir/{wildcards.sample}.hmmsearch.tblout"
        hits_raw="$tmpdir/{wildcards.sample}.hits.raw.tsv"

        marker_db="{params.marker_hmm_db}"
        if [ -d "$marker_db" ]; then
            marker_db_combined="$tmpdir/combined_markers.hmm"
            : > "$marker_db_combined"
            find "$marker_db" -maxdepth 1 -type f \\( -name '*.hmm' -o -name '*.HMM' \\) | sort | while read -r hmm; do
                cat "$hmm" >> "$marker_db_combined"
                printf '\n' >> "$marker_db_combined"
            done
            if [ ! -s "$marker_db_combined" ]; then
                echo "No .hmm/.HMM files found under marker HMM directory: $marker_db" >> {log}
                exit 1
            fi
            marker_db="$marker_db_combined"
        fi

        singularity exec {params.container} hmmsearch \
            --noali \
            --cpu {threads} \
            -E {params.hmmer_evalue} \
            --tblout "$tbl" \
            "$marker_db" {input.faa} >> {log} 2>&1

        # hmmsearch tblout columns used here:
        #   $1 = target gene_id (predicted protein)
        #   $3 = query marker_name (HMM profile)
        #   $5 = full sequence E-value
        awk 'BEGIN{{OFS="\t"}} !/^#/ {{print $1, $3, $5}}' "$tbl" > "$hits_raw"

        if [ ! -s "$hits_raw" ]; then
            echo "No marker hits were generated by hmmscan." >> {log}
            exit 1
        fi

        # Convert protein gene IDs back to contig IDs by removing terminal _[0-9]+.
        # Keep best (lowest e-value) hit per contig_id + marker pair.
        sort -k1,1 -k2,2 -k3,3g "$hits_raw" | \
            awk 'BEGIN{{OFS="\t"}} {{g=$1; sub(/_[0-9]+$/, "", g); key=g"\t"$2; if(!(key in seen)) {{seen[key]=1; print g, $2, $3}}}}' > {output.hmm}
        touch {output.done}
        """

rule magscot_refine_bins:
    """
    Run MAGScoT using bin/contig/set input assembled from metabat2, maxbin2, and concoct bins.
    Marker mapping (--hmm) is generated by magscot_generate_markers.
    """
    input:
        table=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/{{sample}}.input.tsv",
        table_done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/contig_bin_map.done",
        hmm=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/{{sample}}.markers.tsv",
        markers_done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/markers.done"
    output:
        done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/magscot.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot",
        container="containers/magscot.sif",
        out_base=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.sample}/bin_refinement/magscot/" + config.get("magscot", {}).get("out", "MAGScoT"),
        profile=lambda wildcards: config.get("magscot", {}).get("profile", "bac120+ar53"),
        bin_separator=lambda wildcards: config.get("magscot", {}).get("bin_separator", "cleanbin"),
        score_a=lambda wildcards: str(config.get("magscot", {}).get("score_a", 1)),
        score_b=lambda wildcards: str(config.get("magscot", {}).get("score_b", 0.5)),
        score_c=lambda wildcards: str(config.get("magscot", {}).get("score_c", 0.5)),
        max_cont=lambda wildcards: str(config.get("magscot", {}).get("max_cont", 1)),
        threshold=lambda wildcards: str(config.get("magscot", {}).get("threshold", 0.5)),
        min_markers=lambda wildcards: str(config.get("magscot", {}).get("min_markers", 25)),
        min_sharing=lambda wildcards: str(config.get("magscot", {}).get("min_sharing", 0.8)),
        n_iterations=lambda wildcards: str(config.get("magscot", {}).get("n_iterations", 2)),
        score_only=lambda wildcards: bool(config.get("magscot", {}).get("score_only", False)),
        skip_merge_bins=lambda wildcards: bool(config.get("magscot", {}).get("skip_merge_bins", False))
    threads: config["threads"]
    log: "logs/magscot_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        if [ ! -f "{input.hmm}" ]; then
            echo "MAGScoT marker mapping file not found: {input.hmm}" >> {log}
            exit 1
        fi

        if [ ! -s {input.table} ]; then
            echo "MAGScoT contig-to-bin table is empty: {input.table}" >> {log}
            exit 1
        fi

        score_only_flag=""
        skip_merge_flag=""
        profile_arg="{params.profile}"
        if [ "$profile_arg" = "bac120+ar53" ]; then
            profile_arg="default"
        fi
        if [ "{params.score_only}" = "True" ]; then
            score_only_flag="--score_only"
        fi
        if [ "{params.skip_merge_bins}" = "True" ]; then
            skip_merge_flag="--skip_merge_bins"
        fi

        patched_script="{params.outdir}/MAGScoT.patched.R"
        singularity exec {params.container} cp /opt/MAGScoT/MAGScoT.R "$patched_script"

        # Runtime patch for missing funr package in container.
        sed -i 's/suppressMessages(library(funr))/# funr not installed in container; replaced by static script path/' "$patched_script"
        sed -i 's#funr::get_script_path()#"/opt/MAGScoT"#g' "$patched_script"

        singularity exec {params.container} Rscript "$patched_script" \
            -i {input.table} \
            --hmm {input.hmm} \
            -p "$profile_arg" \
            -o {params.out_base} \
            --bin_separator {params.bin_separator} \
            -a {params.score_a} \
            -b {params.score_b} \
            -c {params.score_c} \
            --max_cont {params.max_cont} \
            -t {params.threshold} \
            -m {params.min_markers} \
            -s {params.min_sharing} \
            -n {params.n_iterations} \
            $score_only_flag \
            $skip_merge_flag >> {log} 2>&1

        touch {output.done}
        """

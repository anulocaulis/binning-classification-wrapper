# modules/metawrap_cobinning.smk
# Per-sample, per-assembly MetaWRAP binning/refinement.
# Shared helpers/rules come from modules/readmap_prep.smk.
# As part of metawrap binning (CONCOCT,METABAT2, Maxbin2):

#     Map ALL short reads against each assembly, even long-read only (S1 Flye, S2 Flye, …, S12 IDBAUD)
#     Bin each assembly with those read mappings

# Metawrap bin refinement
# Evaluate binning success with CheckM2

# ---
# Step 2: Bin each sample assembly with all trimmed PolyG-filtered Illumina reads
# ---
rule metawrap_binning_per_assembly:
    input:
        assembly=lambda wildcards: resolve_assembly_path(wildcards.sample, wildcards.assembly_subdir),
        all_r1=expand("trimmed_reads/{sample}/{sample}_R1.fastq.gz", sample=ALL_SAMPLES),
        all_r2=expand("trimmed_reads/{sample}/{sample}_R2.fastq.gz", sample=ALL_SAMPLES)
    output:
        metabat2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/metabat2_bins"),
        maxbin2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/maxbin2_bins"),
        concoct=directory(f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/concoct_bins"),
        maxbin2_status=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/maxbin2.status",
        done=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/metawrap.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}",
        container=METAWRAP_CONTAINER,
        metawrap_exec=config.get("metawrap_executable", "/usr/local/bin/metawrap"),
        min_contig_len=1000,
        retry_min_contig_len=2500
    threads: lambda wildcards: cobinning_threads(wildcards.assembly_subdir)
    resources:
        runtime=lambda wildcards: cobinning_runtime_minutes(wildcards.assembly_subdir),
        mem_mb=lambda wildcards: cobinning_mem_mb(wildcards.assembly_subdir)
    log: "logs/binning_metawrap_{sample}_{assembly_subdir}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        READS_ARGS=""
        i=0
        for r1 in {input.all_r1}; do
            set -- {input.all_r2}
            j=0
            r2=""
            for cand in "$@"; do
                if [ "$j" -eq "$i" ]; then
                    r2="$cand"
                    break
                fi
                j=$((j+1))
            done
            [ -n "$r2" ] || continue
            READS_ARGS="$READS_ARGS $r1 $r2"
            i=$((i+1))
        done

        # MetaWRAP returns non-zero when MaxBin2 fails, even if MetaBAT2/CONCOCT succeeded.
        set +e
        singularity exec {params.container} {params.metawrap_exec} binning \
            -o {params.outdir} \
            -t {threads} \
            -a {input.assembly} \
            --metabat2 --maxbin2 --concoct \
            -l {params.min_contig_len} \
            $READS_ARGS >> {log} 2>&1
        metawrap_rc=$?
        set -e

        if [ ! -d {output.metabat2} ] || [ ! -d {output.concoct} ]; then
            echo "MetaWRAP exited with code $metawrap_rc and required metabat2/concoct output directories are missing." >> {log}
            exit 1
        fi

        if [ "$metawrap_rc" -ne 0 ]; then
            echo "MetaWRAP exited with non-zero status ($metawrap_rc), but required metabat2/concoct outputs exist. Continuing with fallback MaxBin2 handling." >> {log}
        fi

        maxbin2_has_bins=0
        if [ -d {output.maxbin2} ] && ls {output.maxbin2}/*.fa {output.maxbin2}/*.fasta >/dev/null 2>&1; then
            maxbin2_has_bins=1
        fi

        if [ "$maxbin2_has_bins" -eq 1 ]; then
            echo "present" > {output.maxbin2_status}
            rm -f {output.maxbin2}/MAXBIN2_MISSING.txt
            echo "MaxBin2 bins detected for {wildcards.sample}/{wildcards.assembly_subdir}." >> {log}
        else
            retry_outdir="{params.outdir}/retry_maxbin2_min{params.retry_min_contig_len}"
            rm -rf "$retry_outdir"

            echo "MaxBin2 bins missing after first pass; retrying with min_contig_len={params.retry_min_contig_len}." >> {log}
            set +e
            singularity exec {params.container} {params.metawrap_exec} binning \
                -o "$retry_outdir" \
                -t {threads} \
                -a {input.assembly} \
                --metabat2 --maxbin2 --concoct \
                -l {params.retry_min_contig_len} \
                $READS_ARGS >> {log} 2>&1
            retry_rc=$?
            set -e

            retry_maxbin2_has_bins=0
            if [ -d "$retry_outdir/maxbin2_bins" ] && ls "$retry_outdir"/maxbin2_bins/*.fa "$retry_outdir"/maxbin2_bins/*.fasta >/dev/null 2>&1; then
                retry_maxbin2_has_bins=1
            fi

            if [ "$retry_maxbin2_has_bins" -eq 1 ]; then
                rm -rf {output.maxbin2}
                mv "$retry_outdir/maxbin2_bins" {output.maxbin2}
                echo "present" > {output.maxbin2_status}
                rm -f {output.maxbin2}/MAXBIN2_MISSING.txt
                echo "MaxBin2 bins recovered on retry (rc=$retry_rc) with min_contig_len={params.retry_min_contig_len}." >> {log}
            else
                rm -rf {output.maxbin2}
                mkdir -p {output.maxbin2}
                echo "missing" > {output.maxbin2_status}
                printf "MaxBin2 did not produce bins for %s/%s after retries at min_contig_len=%s and %s. Refinement will use metabat2+concoct only.\n" "{wildcards.sample}" "{wildcards.assembly_subdir}" "{params.min_contig_len}" "{params.retry_min_contig_len}" > {output.maxbin2}/MAXBIN2_MISSING.txt
                echo "MaxBin2 bins still missing for {wildcards.sample}/{wildcards.assembly_subdir} after retry (rc=$retry_rc)." >> {log}
            fi

            rm -rf "$retry_outdir"
        fi

        touch {output.done}
        """


# ---
# Step 3: Refine each sample assembly binning outputs
# ---
rule metawrap_refinement_per_assembly:
    input:
        metabat2=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/metabat2_bins",
        maxbin2=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/maxbin2_bins",
        concoct=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/concoct_bins",
        maxbin2_status=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/maxbin2.status",
        metawrap_done=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{{assembly_subdir}}/metawrap.done"
    output:
        refined_bins=directory(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/by_assembly/{{assembly_subdir}}/metawrap_50_10_bins"),
        done=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/by_assembly/{{assembly_subdir}}/refinement.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/by_assembly/{{assembly_subdir}}",
        container=METAWRAP_CONTAINER,
        metawrap_exec=config.get("metawrap_executable", "/usr/local/bin/metawrap"),
        checkm_data_dir=config.get("checkm_db", ""),
        completeness=50,
        contamination=10
    threads: config["threads"]
    log: "logs/binning_refinement_{sample}_{assembly_subdir}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}

        BIN_ARGS="-A {input.metabat2} -B {input.concoct}"
        if [ -f "{input.maxbin2}/MAXBIN2_MISSING.txt" ] || [ "$(cat {input.maxbin2_status})" = "missing" ]; then
            echo "MaxBin2 missing for {wildcards.sample}/{wildcards.assembly_subdir}; running bin_refinement with metabat2+concoct as two bin sets (-A/-B)." >> {log}
        else
            BIN_ARGS="-A {input.metabat2} -B {input.maxbin2} -C {input.concoct}"
            echo "MaxBin2 present for {wildcards.sample}/{wildcards.assembly_subdir}; running bin_refinement with all three bin sets." >> {log}
        fi

        singularity exec {params.container} sh -c "set -euo pipefail; \
            PATH=/opt/conda/envs/metawrap/bin:\$PATH; \
            if [ -n '{params.checkm_data_dir}' ] && [ -d '{params.checkm_data_dir}' ]; then \
                export CHECKM_DATA_PATH='{params.checkm_data_dir}'; \
                export CHECKM_CONFIG_FILE='{params.checkm_data_dir}/DATA_CONFIG'; \
            fi; \
            {params.metawrap_exec} bin_refinement \
                -o {params.outdir} \
                -t {threads} \
                $BIN_ARGS \
                -c {params.completeness} \
                -x {params.contamination}" >> {log} 2>&1

        if [ -d "{params.outdir}/binsO" ]; then
            [ ! -d "{output.refined_bins}" ] && mv "{params.outdir}/binsO" "{output.refined_bins}"
        else
            refined_src=""
            for cand in binsABC binsAB binsAC binsBC binsA binsB binsC; do
                if [ -d "{params.outdir}/work_files/$cand" ]; then
                    refined_src="{params.outdir}/work_files/$cand"
                    break
                fi
            done

            if [ -n "$refined_src" ]; then
                rm -rf "{output.refined_bins}"
                cp -r "$refined_src" "{output.refined_bins}"
            else
                echo "No refined bin directory found (expected binsO or work_files/bins*)." >> {log}
                exit 1
            fi
        fi

        touch {output.done}
        """


# ---
# Step 4: Canonicalize primary assembly outputs for downstream modules
# ---
rule metawrap_select_primary_outputs:
    input:
        metabat2=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{PRIMARY_BINNING_SUBDIR}/metabat2_bins",
        maxbin2=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{PRIMARY_BINNING_SUBDIR}/maxbin2_bins",
        concoct=f"{OUTPUT_DIR}/{{sample}}/binning/by_assembly/{PRIMARY_BINNING_SUBDIR}/concoct_bins",
        refined_bins=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/by_assembly/{PRIMARY_BINNING_SUBDIR}/metawrap_50_10_bins"
    output:
        metabat2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/metabat2_bins"),
        maxbin2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2_bins"),
        concoct=directory(f"{OUTPUT_DIR}/{{sample}}/binning/concoct_bins"),
        refined_bins=directory(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins")
    log: "logs/metawrap_select_primary_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.metabat2}) logs
        : > {log}

        rm -rf {output.metabat2} {output.maxbin2} {output.concoct} {output.refined_bins}

        ln -s {input.metabat2} {output.metabat2}
        ln -s {input.maxbin2} {output.maxbin2}
        ln -s {input.concoct} {output.concoct}
        ln -s {input.refined_bins} {output.refined_bins}

        echo "Linked primary ({PRIMARY_BINNING_SUBDIR}) binning outputs to canonical paths." >> {log}
        """


# ---
# Step 5: CheckM2 bin quality assessment on canonical refined bins
# ---
rule checkm2_quality:
    input:
        bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
    output:
        report=f"{OUTPUT_DIR}/{{sample}}/checkm2/quality_report.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/checkm2",
        container=MAGQUAL_CONTAINER,
        executable="checkm2",
        database=lambda wildcards: config.get("checkm2_db", "")
    threads: config["threads"]
    log: "logs/checkm2_quality_{sample}.log"
    shell:
        """
        if [ -z "{params.database}" ] || [ ! -e "{params.database}" ]; then
            echo "CheckM2 requires config.checkm2_db to point to a downloaded DIAMOND database file." >&2
            exit 1
        fi

        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        singularity exec {params.container} sh -lc 'PATH=/opt/conda/envs/magqual/bin:$PATH; {params.executable} predict \
            --input {input.bins_dir} \
            --output-directory {params.outdir} \
            --threads {threads} \
            --database_path {params.database} \
            --extension .fa' 2>> {log}

        if [ ! -f {output.report} ]; then
            echo "CheckM2 did not produce quality_report.tsv at {output.report}" >&2
            exit 1
        fi
        """

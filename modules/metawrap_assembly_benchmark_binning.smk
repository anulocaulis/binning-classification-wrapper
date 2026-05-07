# This module compares assemblies to report the genome-level outputs
# Make all 7 assemblies for S1 S2 S5(check)
# As part of metawrap binning (CONCOCT,METABAT2, Maxbin2):
# Map all S1 S2 S5 Illumina reads against each assembly that has long reads (S1 Flye, S2 Flye, …, S5 IDBAUD)
# Calculate fraction of reads mapped to each assembly as a measure of assembly success
# Bin each assembly with those read mappings
# Bin refinement module
# Evaluate binning success with CheckM2

# S1, S2, S5 comparison
rule benchmark_binning:
    """
    Bins only the assemblies that have long reads (S1 Flye, S2 Flye, S5 IDBAUD) using metawrap binning
    with the corresponding Illumina reads for each sample. Then runs CheckM2 on the refined bins to evaluate their quality.
    inputs:
        - Assemblies: S1 Flye, S2 Flye, S5 IDBAUD
        - Illumina reads for S1, S2, S5
    outputs:
        - Bins for each assembly and sample combination
        - CheckM2 reports for each set of refined bins
    logs:
        - Logs metawrap binning and CheckM2 runs at logs/cobinning_by_assembly/{{assembly_subdir}}/metawrap.log

    """
    input:
        assembly=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/work/combined.prefixed.min1000.fa",
        staged_r1=lambda wildcards: [
            f"{OUTPUT_DIR}/{sample}/binning/work/{sample}_1.fastq"
            for sample in samples_with_assembly(wildcards.assembly_subdir)
        ],
        staged_r2=lambda wildcards: [
            f"{OUTPUT_DIR}/{sample}/binning/work/{sample}_2.fastq"
            for sample in samples_with_assembly(wildcards.assembly_subdir)
        ]
    output:
        metabat2=directory(f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/metabat2_bins"),
        maxbin2=directory(f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/maxbin2_bins"),
        concoct=directory(f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/concoct_bins"),
        maxbin2_status=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/maxbin2.status",
        done=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/benchmark_metawrap.done"
    params:
        outdir=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap",
        read_pairs=cobinning_staged_read_pairs,
        container=METAWRAP_CONTAINER,
        metawrap_exec=config.get("metawrap_executable", "/usr/local/bin/metawrap"),
        min_contig_len=1000,
        retry_min_contig_len=2500
    threads: lambda wildcards: cobinning_threads(wildcards.assembly_subdir)
    resources:
        runtime=lambda wildcards: cobinning_runtime_minutes(wildcards.assembly_subdir),
        mem_mb=lambda wildcards: cobinning_mem_mb(wildcards.assembly_subdir)
    log:
        f"logs/cobinning_by_assembly/{{assembly_subdir}}/metawrap.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        : > {log}

        READS_ARGS=""
        for pair in {params.read_pairs}; do
            r1="${{pair%%::*}}"
            r2="${{pair#*::}}"
            READS_ARGS="$READS_ARGS $r1 $r2"
        done

        echo "Starting metawrap binning for assembly {input.assembly} at $(date)" >> {log}

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
            echo "MaxBin2 bins detected for {wildcards.assembly_subdir}." >> {log}
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
                printf "MaxBin2 did not produce bins for %s after retries at min_contig_len=%s and %s. Refinement will use metabat2+concoct only.\n" "{wildcards.assembly_subdir}" "{params.min_contig_len}" "{params.retry_min_contig_len}" > {output.maxbin2}/MAXBIN2_MISSING.txt
                echo "MaxBin2 bins still missing for {wildcards.assembly_subdir} after retry (rc=$retry_rc)." >> {log}
            fi

            rm -rf "$retry_outdir"
        fi

        echo "Finished metawrap binning for assembly {input.assembly} at $(date)" >> {log}
        touch {output.done}
        """

rule benchmark_refinement:
    input:
        metabat2=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/metabat2_bins",
        maxbin2=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/maxbin2_bins",
        concoct=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/concoct_bins",
        maxbin2_status=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/maxbin2.status",
        pooled_done=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_metawrap/benchmark_metawrap.done"
    output:
        refined_bins=directory(f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/metawrap_50_10_bins"),
        done=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/benchmark_refinement.done"
    params:
        outdir=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement",
        container=METAWRAP_CONTAINER,
        metawrap_exec=config.get("metawrap_executable", "/usr/local/bin/metawrap"),
        checkm_data_dir=config.get("checkm_db", ""),
        completeness=50,
        contamination=10
    threads: config["threads"]
    log: "logs/cobinning_by_assembly/{assembly_subdir}/benchmark_refinement.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        : > {log}

        BIN_ARGS="-A {input.metabat2} -B {input.concoct}"
        if [ -f "{input.maxbin2}/MAXBIN2_MISSING.txt" ] || [ "$(cat {input.maxbin2_status})" = "missing" ]; then
            echo "MaxBin2 missing for {wildcards.assembly_subdir}; running bin_refinement with metabat2+concoct." >> {log}
        else
            BIN_ARGS="-A {input.metabat2} -B {input.maxbin2} -C {input.concoct}"
            echo "MaxBin2 present for {wildcards.assembly_subdir}; running bin_refinement with all three bin sets." >> {log}
        fi

        singularity exec {params.container} bash -lc "set -euo pipefail; \
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
                echo "No refined benchmark bin directory found (expected binsO or work_files/bins*)." >> {log}
                exit 1
            fi
        fi

        touch {output.done}
        """

rule benchmark_checkm2:
    input:
        bins_dir=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/metawrap_50_10_bins",
        refinement_done=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/benchmark_refinement.done"
    output:
        report=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_checkm2/quality_report.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_checkm2",
        container=MAGQUAL_CONTAINER,
        executable="checkm2",
        database=lambda wildcards: config.get("checkm2_db", "")
    threads: config["threads"]
    log: "logs/cobinning_by_assembly/{assembly_subdir}/benchmark_checkm2.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        : > {log}

        if [ -z "{params.database}" ] || [ ! -e "{params.database}" ]; then
            echo "CheckM2 requires config.checkm2_db to point to a downloaded DIAMOND database file." >> {log}
            exit 1
        fi

        rm -rf {params.outdir}
        mkdir -p {params.outdir} $(dirname {log})

        singularity exec {params.container} sh -lc 'PATH=/opt/conda/envs/magqual/bin:$PATH; {params.executable} predict \
            --input {input.bins_dir} \
            --output-directory {params.outdir} \
            --threads {threads} \
            --database_path {params.database} \
            --extension .fa' >> {log} 2>&1

        if [ ! -f {output.report} ]; then
            echo "CheckM2 did not produce quality_report.tsv at {output.report}" >> {log}
            exit 1
        fi
        """

rule benchmark_select_hqmags:
    input:
        bins_dir=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/metawrap_50_10_bins",
        report=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_checkm2/quality_report.tsv"
    output:
        hqmags=directory(f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/HQMAGs"),
        done=f"{OUTPUT_DIR}/cobinning/by_assembly/{{assembly_subdir}}/benchmark_refinement/hqmags.done"
    log: "logs/cobinning_by_assembly/{assembly_subdir}/benchmark_hqmags.log"
    shell:
        """
        mkdir -p {output.hqmags} $(dirname {log})
        : > {log}
        rm -rf {output.hqmags}
        mkdir -p {output.hqmags}

        awk -F'\t' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    header = tolower($i)
                    gsub(/[^a-z0-9]/, "", header)
                    if (header == "completeness") comp_col = i
                    if (header == "contamination") cont_col = i
                }
                next
            }
            comp_col && cont_col && ($comp_col + 0) > 90 && ($cont_col + 0) < 5 { print $1 }
        ' {input.report} | while IFS= read -r bin_name; do
            [ -n "$bin_name" ] || continue

            src=""
            for cand in \
                "{input.bins_dir}/$bin_name" \
                "{input.bins_dir}/$bin_name.fa" \
                "{input.bins_dir}/$bin_name.fasta"; do
                if [ -f "$cand" ]; then
                    src="$cand"
                    break
                fi
            done

            if [ -n "$src" ]; then
                cp "$src" {output.hqmags}/
                echo "Selected HQMAG: $(basename "$src")" >> {log}
            else
                echo "No source bin found for HQ candidate: $bin_name" >> {log}
            fi
        done

        touch {output.done}
        """
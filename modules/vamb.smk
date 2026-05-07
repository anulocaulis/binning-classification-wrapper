# modules/vamb.smk
# Per-assembly-type pooled VAMB binning across all samples.
# 
# VAMB workflow:
#   1. Concatenate assemblies from all samples for a given assembly_subdir into a single FASTA.
#   2. Map all reads (from all samples) to this pooled catalogue.
#   3. Generate abundance information from BAM files (via TSV).
#   4. Run VAMB binning on the pooled assembly once per assembly_subdir.
#   5. Link outputs back to per-sample paths for compatibility.
#
# Follows VAMB's recommended approach:
#   - Unique contig headers via sample-prefixed naming (e.g., S1C1, S2C2, ...)
#   - Abundance TSV with one column per sample
#   - Sorted BAM files for mapping verification

def vamb_samples_with_assembly(assembly_subdir):
    """
    Return list of samples that have a resolvable assembly for the given assembly_subdir.
    
    **Parameters:**
    - assembly_subdir: Assembly type to query (e.g., 'flye', 'metaspades', 'megahit')
    
    **Returns:**
    - List of sample names (e.g., ['S1', 'S2', ...]) that have valid assemblies for this type.
    
    **Raises:**
    - ValueError if no samples have a resolvable assembly for the given type.
    
    **Logic:**
    Iterates through ALL_SAMPLES and attempts to resolve each sample's assembly path
    for the given assembly_subdir using resolve_assembly_path(). Silently skips samples
    that don't have a valid assembly (FileNotFoundError). Raises an error if the final
    list is empty, preventing silent failures in VAMB rules.
    """
    resolved = []
    for sample in ALL_SAMPLES:
        try:
            resolve_assembly_path(sample, assembly_subdir)
            resolved.append(sample)
        except FileNotFoundError:
            continue
    if not resolved:
        raise ValueError(f"No samples have a resolvable assembly for {assembly_subdir}")
    return resolved


# ---
# Step 1: Concatenate assemblies by assembly_subdir
# ---
rule vamb_concatenate_assemblies:
    """
    Pool assemblies from all samples with a given assembly_subdir.
    Output contig headers follow the format: {Sample}{Contig#} (e.g., S1C1, S1C2, ..., S2C1, ...)
    Minimum contig length is 1000 bp to reduce noise.
    
    **Inputs:**
    - assemblies: One FASTA file per sample for the given assembly_subdir
    
    **Outputs:**
    - fasta: Single concatenated FASTA with all contigs ≥1000 bp, sample-prefixed headers
    
    **Logic:**
    Iterates through each input assembly, numbering samples sequentially (S1, S2, ...).
    For each contig in each assembly, checks if it meets the minimum length threshold before
    output. Rebuilds contig headers to use the sample prefix (S1, S2, etc.) + a sequential 
    counter (C1, C2, ...), ensuring globally unique names across all samples.
    """
    input:
        assemblies=lambda wildcards: [
            resolve_assembly_path(sample, wildcards.assembly_subdir)
            for sample in vamb_samples_with_assembly(wildcards.assembly_subdir)
        ]
    output:
        fasta=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/catalogs/combined.prefixed.min1000.fa"
    params:
        outdir=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/catalogs",
        min_contig_len=1000,
        separator="C"  # Separator between sample name and contig number
    log: "logs/vamb_concatenate_{assembly_subdir}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}
        
        tmpfile=$(mktemp)
        sample_num=0
        
        for asm in {input.assemblies}; do
            sample_num=$((sample_num + 1))
            sample_name="S${{sample_num}}"
            
            if [ ! -f "$asm" ]; then
                echo "Assembly not found: $asm" >> {log}
                rm -f "$tmpfile"
                exit 1
            fi
            
            # Extract FASTA, filter by min_contig_len, and prefix contig headers
            awk -v min={params.min_contig_len} -v sep="{params.separator}" -v sname="$sample_name" '
                /^>/ {{
                    if (len >= min && hdr != "") {{
                        print hdr; print seq
                    }}
                    # Increment counter and build new header (starts at 1)
                    cnt++
                    hdr = ">"sname sep cnt
                    seq = ""
                    len = 0
                }}
                !/^>/ {{
                    seq = seq $0; len += length($0)
                }}
                END {{
                    if (len >= min && hdr != "") {{
                        print hdr; print seq
                    }}
                }}
            ' "$asm" >> "$tmpfile"
        done
        
        mv "$tmpfile" {output.fasta}
        echo "Concatenated $sample_num samples into {output.fasta}." >> {log}
        """


# ---
# Step 2: Map all reads from each sample to the pooled assembly
# ---
rule vamb_map_reads_to_pooled_assembly:
    """
    Map trimmed reads from each sample to the pooled assembly catalogue.
    Outputs a sorted, indexed BAM file.
    
    **Inputs:**
    - fasta: Pooled assembly FASTA (from vamb_concatenate_assemblies)
    - reads: Interleaved trimmed FASTQ for this sample
    
    **Outputs:**
    - bam: Coordinate-sorted BAM file
    - bai: BAM index file
    
    **Logic:**
    Uses minimap2 with short-read parameters (-ax sr) to map reads to the catalogue,
    pipes directly to samtools sort (avoiding intermediate SAM), then indexes the result.
    The --split-prefix ensures temporary index files go to /tmp/ rather than output directory.
    Validates that the BAM is non-empty before completing.
    """
    input:
        fasta=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/catalogs/combined.prefixed.min1000.fa",
        reads=f"trimmed_reads/{{sample}}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        bam=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/mapping/{{sample}}.sorted.bam.bai"
    params:
        outdir=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/mapping",
        container=QC_CONTAINER
    threads: config.get("vamb_map_threads", config.get("threads", 8))
    resources:
        mem_mb=lambda wildcards: 8000,
        runtime=lambda wildcards: 120
    log: "logs/vamb_map_{assembly_subdir}_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}
        
        singularity exec {params.container} sh -c "set -euo pipefail; \
            export PATH=/opt/conda/envs/qc_env/bin:/opt/conda/envs/base_tools/bin:/opt/conda/bin:\$PATH; \
            minimap2 -t {threads} -ax sr --split-prefix /tmp/mm2.tmp {input.fasta} {input.reads} 2>> {log} | \
            samtools sort -@ {threads} -m 4G -o {output.bam} 2>> {log}; \
            samtools index -@ {threads} {output.bam} 2>> {log}" 2>&1 | tee -a {log}
        
        if [ ! -f "{output.bam}" ] || [ ! -s "{output.bam}" ]; then
            echo "Error: BAM file not created or is empty." >> {log}
            exit 1
        fi
        """


# ---
# Step 3: Generate abundance TSV from BAM files
# ---
rule vamb_build_abundance_tsv:
    """
    Build abundance TSV with contig names and per-sample abundance columns.
    Abundance is calculated as: (mapped reads) / (reference length)
    
    **Inputs:**
    - fasta: Pooled assembly FASTA
    - bams: One sorted BAM per sample
    - bais: Corresponding BAM indices
    
    **Outputs:**
    - tsv: Tab-separated file with header (contigname + sample columns) and per-contig abundances
    
    **Logic:**
    For each sample, samtools idxstats reports per-reference statistics (length, mapped, unmapped).
    Abundance per contig is calculated as mapped_count / reference_length.
    Each sample's abundance is joined into an accumulated table using the 'join' command with:
    - '-a 1' to keep all contigs from the left (accumulated) table, even if missing in current sample
    - '-e 0' to fill missing abundances as 0
    - This iterative approach builds a matrix where rows are contigs and columns are samples.
    """
    input:
        fasta=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/catalogs/combined.prefixed.min1000.fa",
        bams=lambda wildcards: [
            f"{OUTPUT_DIR}/vamb/by_assembly/{wildcards.assembly_subdir}/mapping/{sample}.sorted.bam"
            for sample in vamb_samples_with_assembly(wildcards.assembly_subdir)
        ],
        bais=lambda wildcards: [
            f"{OUTPUT_DIR}/vamb/by_assembly/{wildcards.assembly_subdir}/mapping/{sample}.sorted.bam.bai"
            for sample in vamb_samples_with_assembly(wildcards.assembly_subdir)
        ]
    output:
        tsv=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/combined.abundance.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}",
        container=QC_CONTAINER,
        samples=lambda wildcards: " ".join(vamb_samples_with_assembly(wildcards.assembly_subdir))
    threads: config.get("vamb_abundance_threads", 4)
    log: "logs/vamb_abundance_{assembly_subdir}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        : > {log}
        
        tmpdir=$(mktemp -d)
        
        # Extract all contig names from FASTA
        awk '/^>/ {{h=substr($0,2); sub(/[ \t].*$/, "", h); print h}}' {input.fasta} | sort -u > "$tmpdir/contigs.tsv"
        cp "$tmpdir/contigs.tsv" "$tmpdir/table.tsv"
        
        # For each sample, compute per-contig abundance and join
        sample_count=0
        for s in {params.samples}; do
            sample_count=$((sample_count+1))
            bam="{OUTPUT_DIR}/vamb/by_assembly/{wildcards.assembly_subdir}/mapping/${{s}}.sorted.bam"
            
            if [ ! -f "$bam" ]; then
                echo "BAM file not found: $bam" >> {log}
                rm -rf "$tmpdir"
                exit 1
            fi
            
            # Use samtools idxstats to get per-reference abundance
            singularity exec {params.container} samtools idxstats "$bam" 2>> {log} | \
                awk 'BEGIN{{OFS="\t"}} $1 != "*" {{
                    ref_len=$2+0
                    mapped=$3+0
                    unmapped=$4+0
                    abund=(ref_len>0 ? mapped/ref_len : 0)
                    print $1, abund
                }}' | \
                sort -k1,1 > "$tmpdir/${{s}}.tsv"
            
            # Join with accumulated abundance table
            join -t $'\\t' -a 1 -e 0 "$tmpdir/table.tsv" "$tmpdir/${{s}}.tsv" > "$tmpdir/table.next.tsv"
            mv "$tmpdir/table.next.tsv" "$tmpdir/table.tsv"
        done
        
        # Write final TSV with header
        printf 'contigname\\t{params.samples}\n' | sed 's/ /\\t/g' > {output.tsv}
        cat "$tmpdir/table.tsv" >> {output.tsv}
        
        contig_count=$(tail -n +2 {output.tsv} | wc -l)
        if [ "$contig_count" -lt 1 ]; then
            echo "Error: Abundance TSV has no contigs (only header)." >> {log}
            rm -rf "$tmpdir"
            exit 1
        fi
        
        echo "Generated abundance TSV with $contig_count contigs and $sample_count samples." >> {log}
        rm -rf "$tmpdir"
        """


# ---
# Step 4: Run VAMB binning on pooled assemblies
# ---
rule vamb_binning:
    """
    Run VAMB clustering on the pooled assembly using abundance TSV.
    Produces bins (clusters) stored in outdir/clusters.tsv.
    
    **Inputs:**
    - fasta: Pooled assembly FASTA
    - abundance: Abundance TSV with per-sample columns
    
    **Outputs:**
    - bins: Output directory containing all VAMB results
    - clusters_tsv: Primary output file listing contig-to-cluster assignments
    - done: Marker file indicating successful completion
    
    **Logic:**
    Invokes 'vamb bin default' with the provided FASTA and abundance matrix.
    The --norefcheck flag may be set via config to skip reference consistency checks.
    Validates that clusters.tsv exists and is non-empty before completing, ensuring
    successful clustering occurred. Removes any pre-existing output directory to avoid
    contamination from prior runs.
    """
    input:
        fasta=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/catalogs/combined.prefixed.min1000.fa",
        abundance=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/combined.abundance.tsv"
    output:
        bins=directory(f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/vamb_bins"),
        clusters_tsv=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/vamb_bins/clusters.tsv",
        done=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/vamb_binning.done"
    params:
        outdir=f"{OUTPUT_DIR}/vamb/by_assembly/{{assembly_subdir}}/vamb_bins",
        container=VAMB_CONTAINER,
        min_contig_len=lambda wildcards: config.get("vamb", {}).get("min_contig_len", 2000),
        threads=lambda wildcards: config.get("vamb", {}).get("threads", config.get("threads", 8)),
        norefcheck=lambda wildcards: config.get("vamb", {}).get("norefcheck", False)
    threads: config.get("vamb_threads", config.get("threads", 8))
    resources:
        mem_mb=lambda wildcards: config.get("vamb", {}).get("mem_mb", 32000),
        runtime=lambda wildcards: config.get("vamb", {}).get("runtime", 360)
    log: "logs/vamb_binning_{assembly_subdir}.log"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir} logs
        : > {log}
        
        norefcheck_flag=""
        if [ "{params.norefcheck}" = "True" ]; then
            norefcheck_flag="--norefcheck"
        fi
        
        singularity exec {params.container} vamb bin default \
            --outdir {params.outdir} \
            --fasta {input.fasta} \
            --abundance_tsv {input.abundance} \
            -m {params.min_contig_len} \
            -p {params.threads} \
            $norefcheck_flag 2>&1 | tee {log}
        
        # Verify output
        if [ ! -f "{output.clusters_tsv}" ] || [ ! -s "{output.clusters_tsv}" ]; then
            echo "Error: VAMB finished but clusters.tsv is missing or empty in {params.outdir}" >> {log}
            exit 1
        fi
        
        touch {output.done}
        """


# ---
# Step 5: Link per-assembly VAMB outputs back to per-sample compatibility paths
# ---
rule vamb_select_primary_outputs:
    """
    Link pooled VAMB outputs from PRIMARY_BINNING_SUBDIR to per-sample canonical paths.
    This maintains compatibility with downstream rules that expect per-sample outputs.
    
    **Inputs:**
    - bins: Pooled VAMB bin directory from primary assembly type
    - done: Marker indicating successful VAMB binning completion
    
    **Outputs:**
    - bins: Per-sample symbolic link to the shared pooled bins directory
    - done: Per-sample marker file indicating this sample's outputs are ready
    
    **Logic:**
    Creates a symbolic link from per-sample path to the shared pooled result directory
    (rather than duplicating large bin files). This allows downstream rules to access
    per-sample-styled paths while only storing one actual copy of the bins. The realpath
    call ensures the symlink is absolute, avoiding broken relative paths if directories
    are moved.
    """
    input:
        bins=f"{OUTPUT_DIR}/vamb/by_assembly/{PRIMARY_BINNING_SUBDIR}/vamb_bins",
        done=f"{OUTPUT_DIR}/vamb/by_assembly/{PRIMARY_BINNING_SUBDIR}/vamb_binning.done"
    output:
        bins=directory(f"{OUTPUT_DIR}/{{sample}}/binning/vamb_bins"),
        done=f"{OUTPUT_DIR}/{{sample}}/binning/vamb.done"
    log: "logs/vamb_select_primary_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.bins}) logs
        : > {log}
        
        rm -rf {output.bins}
        ln -s $(realpath {input.bins}) {output.bins}
        echo "Linked primary VAMB bins from {input.bins} to {output.bins}" >> {log}
        
        touch {output.done}
        """

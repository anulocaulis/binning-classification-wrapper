#!/usr/bin/env bash
#SBATCH --job-name=read_fraction_summary
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/read_fraction_summary_%j.log
#SBATCH --error=logs/read_fraction_summary_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mike.beitner@ucdenver.edu
set -euo pipefail

# Summarize read fractions by mapping each sample's reads to each assembler output:
#   /storage/biology/projects/miller-lowry/beitner/data/assemblies/<assembler>
#
# Usage:
#   sbatch summarize_read_fractions.sh [output_tsv]
#   bash   summarize_read_fractions.sh [output_tsv]
#
# Optional tuning:
#   PARALLEL_JOBS=4 MAP_THREADS=4 sbatch summarize_read_fractions.sh

OUTPUT_TSV="${1:-binning_outputs/summary/read_fraction_summary.tsv}"
OUTPUT_DIR="$(dirname "$OUTPUT_TSV")"

ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-/storage/biology/projects/miller-lowry/beitner/data/assemblies}"
DATA_ROOT="${DATA_ROOT:-/storage/biology/projects/miller-lowry/beitner/data}"
CONTAINER="${CONTAINER:-containers/multi_qc.sif}"
TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-16}"
PARALLEL_JOBS="${PARALLEL_JOBS:-4}"
MAP_THREADS="${MAP_THREADS:-$(( TOTAL_CPUS / PARALLEL_JOBS ))}"
MAP_CACHE_DIR="${MAP_CACHE_DIR:-binning_outputs/summary/mapping_by_assembly}"

ASSEMBLERS=(flye idbaud megahit metaconnet metamdbg metaspades metaspades_hybrid)

if [[ "$PARALLEL_JOBS" -lt 1 ]]; then
    PARALLEL_JOBS=1
fi
if [[ "$MAP_THREADS" -lt 1 ]]; then
    MAP_THREADS=1
fi

mkdir -p "$OUTPUT_DIR" "$MAP_CACHE_DIR" logs

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

resolve_assembly_fasta() {
    local sample="$1"
    local assembler="$2"
    local sample_lower
    sample_lower="$(echo "$sample" | tr '[:upper:]' '[:lower:]')"
    local base="${ASSEMBLY_ROOT}/${assembler}"
    local cand

    for cand in \
        "${base}/${sample}.assembly.fasta" \
        "${base}/${sample_lower}.assembly.fasta" \
        "${base}/${sample}.assembly.fa" \
        "${base}/${sample_lower}.assembly.fa" \
        "${base}/${sample}.${assembler}.fasta" \
        "${base}/${sample_lower}.${assembler}.fasta" \
        "${base}/${sample}.${assembler}.fa" \
        "${base}/${sample_lower}.${assembler}.fa"
    do
        if [[ -s "$cand" ]]; then
            echo "$cand"
            return 0
        fi
    done
    return 1
}

find_bins_root() {
    local sample="$1"
    local assembly_subdir="$2"
    local candidate

    for candidate in \
        "binning_outputs/${sample}/bin_refinement/by_assembly/${assembly_subdir}" \
        "binning_outputs/${sample}/bin_refinement/by_assembly/${assembly_subdir}/metawrap_50_10_bins" \
        "binning_outputs/${sample}/bin_refinement/metawrap_50_10_bins"
    do
        if [[ -d "$candidate" ]]; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

process_sample_assembly() {
    local sample="$1"
    local reads="$2"
    local assembly="$3"
    local rows_dir="$4"
    local assembly_subdir="assembly.${assembly}"
    local rowfile="${rows_dir}/${sample}__${assembly}.tsv"
    local task_log="logs/read_fraction_map_${sample}_${assembly}.log"
    local notes="ok"
    local assembly_fasta
    local map_subdir bam bai flagstat idxstats
    local total_reads mapped_reads assembled_fraction
    local binned_mapped_reads binned_fraction_total binned_fraction_mapped
    local bins_root key bins_contigs mapped_by_contig joined

    if ! assembly_fasta="$(resolve_assembly_fasta "$sample" "$assembly")"; then
        printf "%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_assembly\n" "$sample" "$assembly" "$assembly_subdir" > "$rowfile"
        return 0
    fi

    map_subdir="${MAP_CACHE_DIR}/${assembly}"
    mkdir -p "$map_subdir"
    bam="${map_subdir}/${sample}.sorted.bam"
    bai="${bam}.bai"
    flagstat="${map_subdir}/${sample}.flagstat.txt"
    idxstats="${map_subdir}/${sample}.idxstats.txt"

    if [[ ! -s "$bam" || ! -s "$bai" || ! -s "$flagstat" || ! -s "$idxstats" ]]; then
        singularity exec "$CONTAINER" sh -c "set -euo pipefail; \
            export PATH=/opt/conda/envs/qc_env/bin:/opt/conda/envs/base_tools/bin:/opt/conda/bin:\$PATH; \
            minimap2 -t ${MAP_THREADS} -ax sr --split-prefix '${bam}.mmi' '${assembly_fasta}' '${reads}' | \
            samtools sort -@ ${MAP_THREADS} -o '${bam}'; \
            samtools index -@ ${MAP_THREADS} '${bam}'; \
            samtools flagstat -@ ${MAP_THREADS} '${bam}' > '${flagstat}'; \
            samtools idxstats '${bam}' > '${idxstats}'" > "$task_log" 2>&1
    fi

    total_reads="$(awk '/ in total /{print $1; exit}' "$flagstat")"
    mapped_reads="$(awk '/ mapped \(/ && $0 !~ /primary mapped \(/ {print $1; exit}' "$flagstat")"
    if [[ -z "$total_reads" || -z "$mapped_reads" ]]; then
        printf "%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tflagstat_parse_error\n" "$sample" "$assembly" "$assembly_subdir" > "$rowfile"
        return 0
    fi

    assembled_fraction="$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{if(t>0) printf "%.6f", m/t; else print "0.000000"}')"
    binned_mapped_reads=0

    if bins_root="$(find_bins_root "$sample" "$assembly_subdir")"; then
        key="${sample}.${assembly}"
        bins_contigs="$tmpdir/${key}.bins_contigs.tsv"
        mapped_by_contig="$tmpdir/${key}.mapped_by_contig.tsv"
        joined="$tmpdir/${key}.joined.tsv"

        find "$bins_root" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) -print0 \
            | xargs -0 -r awk '/^>/{h=substr($0,2); sub(/[ \t].*$/, "", h); print h}' \
            | sort -u > "$bins_contigs"

        if [[ -s "$bins_contigs" ]]; then
            awk 'BEGIN{OFS="\t"} $1!="*"{print $1, $3}' "$idxstats" | sort -k1,1 > "$mapped_by_contig"
            join -t $'\t' -1 1 -2 1 "$bins_contigs" "$mapped_by_contig" > "$joined" || true
            if [[ -s "$joined" ]]; then
                binned_mapped_reads="$(awk '{s+=$2} END{print s+0}' "$joined")"
            else
                notes="no_binned_contig_hits"
            fi
        else
            notes="no_bin_contigs"
        fi
    else
        notes="missing_bins_dir"
    fi

    binned_fraction_total="$(awk -v b="$binned_mapped_reads" -v t="$total_reads" 'BEGIN{if(t>0) printf "%.6f", b/t; else print "0.000000"}')"
    binned_fraction_mapped="$(awk -v b="$binned_mapped_reads" -v m="$mapped_reads" 'BEGIN{if(m>0) printf "%.6f", b/m; else print "0.000000"}')"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$assembly" "$assembly_subdir" "$total_reads" "$mapped_reads" "$assembled_fraction" "$binned_mapped_reads" "$binned_fraction_total" "$binned_fraction_mapped" "$notes" \
        > "$rowfile"
}

rows_dir="$tmpdir/rows"
mkdir -p "$rows_dir"

active_jobs=0
shopt -s nullglob
for reads in "${DATA_ROOT}"/S*/short_reads.fastq.gz; do
    [[ -s "$reads" ]] || continue
    sample="$(basename "$(dirname "$reads")")"

    for assembly in "${ASSEMBLERS[@]}"; do
        process_sample_assembly "$sample" "$reads" "$assembly" "$rows_dir" &
        active_jobs=$((active_jobs + 1))
        if [[ "$active_jobs" -ge "$PARALLEL_JOBS" ]]; then
            wait -n
            active_jobs=$((active_jobs - 1))
        fi
    done
done
wait

printf "sample\tassembly\tassembly_subdir\ttotal_reads\tmapped_reads\tassembled_fraction\tbinned_mapped_reads\tbinned_fraction_of_total\tbinned_fraction_of_mapped\tnotes\n" > "$OUTPUT_TSV"
cat "$rows_dir"/*.tsv | sort -t $'\t' -k1,1V -k2,2 >> "$OUTPUT_TSV"

grand_total=0
grand_mapped=0
grand_binned=0
rows_count=0
declare -A assembly_totals=()
declare -A assembly_mapped=()
declare -A assembly_binned=()
declare -A assembly_rows=()

while IFS=$'\t' read -r sample assembly assembly_subdir total_reads mapped_reads _assembled_fraction binned_mapped_reads _binned_fraction_total _binned_fraction_mapped _notes; do
    if [[ "$sample" == "sample" || "$sample" == "TOTAL" ]]; then
        continue
    fi
    if [[ "$total_reads" == "NA" || "$mapped_reads" == "NA" || "$binned_mapped_reads" == "NA" ]]; then
        continue
    fi

    grand_total=$((grand_total + total_reads))
    grand_mapped=$((grand_mapped + mapped_reads))
    grand_binned=$((grand_binned + binned_mapped_reads))
    rows_count=$((rows_count + 1))

    assembly_totals["$assembly"]=$(( ${assembly_totals["$assembly"]:-0} + total_reads ))
    assembly_mapped["$assembly"]=$(( ${assembly_mapped["$assembly"]:-0} + mapped_reads ))
    assembly_binned["$assembly"]=$(( ${assembly_binned["$assembly"]:-0} + binned_mapped_reads ))
    assembly_rows["$assembly"]=$(( ${assembly_rows["$assembly"]:-0} + 1 ))
done < "$OUTPUT_TSV"

if [[ "$rows_count" -gt 0 ]]; then
    while IFS= read -r assembly; do
        total_reads="${assembly_totals[$assembly]}"
        mapped_reads="${assembly_mapped[$assembly]}"
        binned_mapped_reads="${assembly_binned[$assembly]}"
        assembly_rows_count="${assembly_rows[$assembly]}"
        assembly_subdir="assembly.${assembly}"

        assembled_fraction="$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{if(t>0) printf "%.6f", m/t; else print "0.000000"}')"
        binned_fraction_total="$(awk -v b="$binned_mapped_reads" -v t="$total_reads" 'BEGIN{if(t>0) printf "%.6f", b/t; else print "0.000000"}')"
        binned_fraction_mapped="$(awk -v b="$binned_mapped_reads" -v m="$mapped_reads" 'BEGIN{if(m>0) printf "%.6f", b/m; else print "0.000000"}')"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "TOTAL" "$assembly" "$assembly_subdir" "$total_reads" "$mapped_reads" "$assembled_fraction" "$binned_mapped_reads" "$binned_fraction_total" "$binned_fraction_mapped" "aggregated_over_${assembly_rows_count}_rows" \
            >> "$OUTPUT_TSV"
    done < <(printf '%s\n' "${!assembly_totals[@]}" | sort)

    total_assembled_fraction="$(awk -v m="$grand_mapped" -v t="$grand_total" 'BEGIN{if(t>0) printf "%.6f", m/t; else print "0.000000"}')"
    total_binned_fraction_total="$(awk -v b="$grand_binned" -v t="$grand_total" 'BEGIN{if(t>0) printf "%.6f", b/t; else print "0.000000"}')"
    total_binned_fraction_mapped="$(awk -v b="$grand_binned" -v m="$grand_mapped" 'BEGIN{if(m>0) printf "%.6f", b/m; else print "0.000000"}')"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "TOTAL" "ALL" "ALL" "$grand_total" "$grand_mapped" "$total_assembled_fraction" "$grand_binned" "$total_binned_fraction_total" "$total_binned_fraction_mapped" "aggregated_over_${rows_count}_rows" \
        >> "$OUTPUT_TSV"
fi

echo "Wrote summary: $OUTPUT_TSV"
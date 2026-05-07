#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <assembly_subdir> <metabat2|concoct>"
    exit 1
fi

assembly_subdir="$1"
binner="$2"

if [[ "$binner" != "metabat2" && "$binner" != "concoct" ]]; then
    echo "Binner must be metabat2 or concoct"
    exit 1
fi

project_root="/storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper"
assemblies_root="/storage/biology/projects/miller-lowry/beitner/data/assemblies"
cd "$project_root"

workdir="binning_outputs/cobinning/by_assembly/${assembly_subdir}/work"
contigs="${workdir}/combined.prefixed.min1000.fa"
outdir="binning_outputs/cobinning/by_assembly/${assembly_subdir}/metawrap"
logfile="logs/manual_${assembly_subdir}_${binner}.log"

declare -a samples
assembly_method=""
case "$assembly_subdir" in
    assembly.idbaud)
        assembly_method="idbaud"
        samples=(S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12)
        ;;
    assembly.megahit)
        assembly_method="megahit"
        samples=(S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12)
        ;;
    assembly.metaspades)
        assembly_method="metaspades"
        samples=(S1 S2 S3 S5 S6 S7 S8 S9 S10 S11 S12)
        ;;
    *)
        echo "Unsupported assembly_subdir: $assembly_subdir"
        exit 1
        ;;
esac

mkdir -p "$workdir" "$outdir" logs

tmp_contigs="${contigs}.tmp_$(date +%Y%m%d_%H%M%S)_$$"
rm -f "$tmp_contigs"
touch "$tmp_contigs"

for s in "${samples[@]}"; do
    s_lower=$(printf "%s" "$s" | tr '[:upper:]' '[:lower:]')
    src_asm="${assemblies_root}/${assembly_method}/${s_lower}.assembly.fasta"
    if [[ ! -s "$src_asm" ]]; then
        echo "Missing assembly for ${s}: $src_asm"
        exit 1
    fi

    awk -v prefix="${s}|" -v min_len="1000" '
        BEGIN { RS=">"; ORS="" }
        NR == 1 { next }
        {
            n = split($0, lines, "\n")
            hdr = lines[1]
            seq = ""
            for (i = 2; i <= n; i++) {
                gsub(/[[:space:]]/, "", lines[i])
                seq = seq lines[i]
            }
            if (length(seq) >= min_len) {
                print ">" prefix hdr "\n"
                for (j = 1; j <= length(seq); j += 80) {
                    print substr(seq, j, 80) "\n"
                }
            }
        }
    ' "$src_asm" >> "$tmp_contigs"
done

if [[ ! -s "$tmp_contigs" ]]; then
    echo "No >=1000bp contigs were produced for ${assembly_subdir} from ${assemblies_root}/${assembly_method}"
    exit 1
fi

mv "$tmp_contigs" "$contigs"

declare -a reads_args
for s in "${samples[@]}"; do
    r1="binning_outputs/${s}/binning/work/${s}_1.fastq"
    r2="binning_outputs/${s}/binning/work/${s}_2.fastq"
    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "Missing reads for ${s}: $r1 or $r2"
        exit 1
    fi
    reads_args+=("$r1" "$r2")
done

tmpout="${outdir}/tmp_${binner}_$(date +%Y%m%d_%H%M%S)_$$"

threads="${SLURM_CPUS_PER_TASK:-32}"

echo "[$(date)] Starting ${assembly_subdir} ${binner} with ${threads} threads" | tee -a "$logfile"

set +e
singularity exec containers/metawrap.sif /usr/local/bin/metawrap binning \
    -o "$tmpout" \
    -t "$threads" \
    -a "$contigs" \
    "--${binner}" \
    -l 1000 \
    "${reads_args[@]}" >> "$logfile" 2>&1
rc=$?
set -e

if [[ $rc -ne 0 ]]; then
    echo "[$(date)] metawrap returned non-zero ($rc) for ${assembly_subdir} ${binner}" | tee -a "$logfile"
fi

produced_dir="${tmpout}/${binner}_bins"
if [[ ! -d "$produced_dir" ]]; then
    echo "[$(date)] Expected output directory missing: $produced_dir" | tee -a "$logfile"
    exit 1
fi

produced_count=$(find "$produced_dir" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | wc -l)
if [[ "$produced_count" -eq 0 ]]; then
    echo "[$(date)] No bins produced in $produced_dir" | tee -a "$logfile"
    exit 1
fi

target_dir="${outdir}/${binner}_bins"
if [[ -d "$target_dir" ]]; then
    mv "$target_dir" "${target_dir}.bak_$(date +%Y%m%d_%H%M%S)"
fi
mv "$produced_dir" "$target_dir"
rm -rf "$tmpout"

metabat2_count=0
concoct_count=0

if [[ -d "${outdir}/metabat2_bins" ]]; then
    metabat2_count=$(find "${outdir}/metabat2_bins" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | wc -l)
fi
if [[ -d "${outdir}/concoct_bins" ]]; then
    concoct_count=$(find "${outdir}/concoct_bins" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | wc -l)
fi

if [[ "$metabat2_count" -gt 0 && "$concoct_count" -gt 0 ]]; then
    mkdir -p "${outdir}/maxbin2_bins"
    echo "missing" > "${outdir}/maxbin2.status"
    printf "MaxBin2 was not run in manual recovery for %s. Refinement will use metabat2+concoct only.\n" "$assembly_subdir" > "${outdir}/maxbin2_bins/MAXBIN2_MISSING.txt"
    touch "${outdir}/cobinning_metawrap.done"
    echo "[$(date)] Finalized ${assembly_subdir}: metabat2=${metabat2_count}, concoct=${concoct_count}" | tee -a "$logfile"
else
    echo "[$(date)] Partial completion for ${assembly_subdir}: metabat2=${metabat2_count}, concoct=${concoct_count}" | tee -a "$logfile"
fi

echo "[$(date)] Completed ${assembly_subdir} ${binner}" | tee -a "$logfile"

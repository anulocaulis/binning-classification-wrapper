#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <assembly_subdir>"
    exit 1
fi

assembly_subdir="$1"
project_root="/storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper"
cd "$project_root"

contigs="binning_outputs/cobinning/by_assembly/${assembly_subdir}/work/combined.prefixed.min1000.fa"
vamb_root="binning_outputs/cobinning/by_assembly/${assembly_subdir}/vamb"
mapping_dir="${vamb_root}/mapping"
outdir="${vamb_root}/vamb_bins"
abundance_tsv="${vamb_root}/combined.abundance.tsv"
logfile="logs/manual_vamb_${assembly_subdir}.log"

qc_container="${QC_CONTAINER_OVERRIDE:-containers/qc_tools_miniconda.sif}"
vamb_container="${VAMB_CONTAINER_OVERRIDE:-containers/vamb.sif}"
threads="${SLURM_CPUS_PER_TASK:-32}"
min_len="${VAMB_MIN_LEN_OVERRIDE:-1000}"

if [[ ! -s "$contigs" ]]; then
    echo "Missing min1000 pooled contigs: $contigs"
    exit 1
fi
if [[ ! -f "$qc_container" ]]; then
    echo "QC container not found: $qc_container"
    exit 1
fi
if [[ ! -f "$vamb_container" ]]; then
    echo "VAMB container not found: $vamb_container"
    exit 1
fi

mkdir -p "$mapping_dir" logs
: > "$logfile"

echo "[$(date)] Starting manual VAMB for ${assembly_subdir} (threads=${threads}, min_len=${min_len})" | tee -a "$logfile"

mapfile -t samples < <(
    awk '/^>/{h=substr($0,2); split(h,a,"|"); if (a[1] ~ /^S[0-9]+$/) print a[1]}' "$contigs" | sort -Vu
)

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "No sample prefixes (S1, S2, ...) were detected in $contigs" | tee -a "$logfile"
    exit 1
fi

echo "[$(date)] Samples detected in contigs: ${samples[*]}" | tee -a "$logfile"

for s in "${samples[@]}"; do
    r1="binning_outputs/${s}/binning/work/${s}_1.fastq"
    r2="binning_outputs/${s}/binning/work/${s}_2.fastq"
    bam="${mapping_dir}/${s}.sorted.bam"
    bai="${bam}.bai"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "Missing reads for ${s}: $r1 or $r2" | tee -a "$logfile"
        exit 1
    fi

    if [[ -s "$bam" && -s "$bai" ]]; then
        echo "[$(date)] Reusing existing mapping for ${s}" >> "$logfile"
        continue
    fi

    echo "[$(date)] Mapping ${s} to pooled contigs" | tee -a "$logfile"
    singularity exec "$qc_container" bash -lc "set -euo pipefail; \
        export PATH=/opt/conda/envs/qc_env/bin:/opt/conda/envs/base_tools/bin:/opt/conda/bin:\$PATH; \
        minimap2 -t ${threads} -ax sr --split-prefix ${bam}.mmi '${contigs}' '${r1}' '${r2}' | \
        samtools sort -@ ${threads} -o '${bam}'; \
        samtools index -@ ${threads} '${bam}'" >> "$logfile" 2>&1

    if [[ ! -s "$bam" || ! -s "$bai" ]]; then
        echo "Mapping/index failed for ${s}" | tee -a "$logfile"
        exit 1
    fi
done

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

awk '/^>/{h=substr($0,2); sub(/[ \t].*$/, "", h); print h}' "$contigs" | sort -u > "${tmpdir}/contigs.tsv"
cp "${tmpdir}/contigs.tsv" "${tmpdir}/table.tsv"

for s in "${samples[@]}"; do
    bam="${mapping_dir}/${s}.sorted.bam"
    singularity exec "$qc_container" samtools idxstats "$bam" | \
        awk 'BEGIN{OFS="\t"} $1 != "*" {len=$2+0; mapped=$3+0; abund=(len>0 ? mapped/len : 0); print $1, abund}' | \
        sort -k1,1 > "${tmpdir}/${s}.tsv"

    join -t $'\t' -a 1 -e 0 "${tmpdir}/table.tsv" "${tmpdir}/${s}.tsv" > "${tmpdir}/table.next.tsv"
    mv "${tmpdir}/table.next.tsv" "${tmpdir}/table.tsv"
done

printf 'contigname\t%s\n' "$(printf '%s\t' "${samples[@]}" | sed 's/\t$//')" > "$abundance_tsv"
cat "${tmpdir}/table.tsv" >> "$abundance_tsv"

if [[ $(wc -l < "$abundance_tsv") -le 1 ]]; then
    echo "Combined abundance TSV has no contig rows: $abundance_tsv" | tee -a "$logfile"
    exit 1
fi

tmpout="${vamb_root}/vamb_bins.tmp_$(date +%Y%m%d_%H%M%S)_$$"
rm -rf "$tmpout"

echo "[$(date)] Running VAMB bin default" | tee -a "$logfile"
singularity exec "$vamb_container" vamb bin default \
    --outdir "$tmpout" \
    --fasta "$contigs" \
    --abundance_tsv "$abundance_tsv" \
    -m "$min_len" \
    -p "$threads" >> "$logfile" 2>&1

if [[ ! -s "${tmpout}/clusters.tsv" ]]; then
    echo "VAMB finished but clusters.tsv is missing in ${tmpout}" | tee -a "$logfile"
    exit 1
fi

if [[ -d "$outdir" ]]; then
    mv "$outdir" "${outdir}.bak_$(date +%Y%m%d_%H%M%S)"
fi
mv "$tmpout" "$outdir"

touch "${vamb_root}/cobinning_vamb.done"

echo "[$(date)] Completed manual VAMB for ${assembly_subdir}" | tee -a "$logfile"

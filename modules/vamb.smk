# modules/vamb.smk
# VAMB binning from assembly contigs + abundance table (no taxonomic annotation)

def get_vamb_assembly(wildcards):
    return config["input_reads"]["assembly"].format(sample=wildcards.sample)


rule vamb_generate_abundance_tsv:
    """
    Build a VAMB-compatible abundance TSV from SR mapped BAM depth.
    Output format header: contigname <tab> sample
    """
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/sr/mapping/short_reads.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/sr/mapping/short_reads.sorted.bam.bai"
    output:
        tsv=f"{OUTPUT_DIR}/{{sample}}/binning/vamb/{{sample}}.abundance.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/binning/vamb",
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/vamb_abundance_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        printf 'contigname\t{wildcards.sample}\n' > {output.tsv}
        singularity exec {params.container} samtools idxstats {input.bam} | \
            awk 'BEGIN{{OFS="\t"}} $1 != "*" {{len=$2+0; mapped=$3+0; abund=(len>0 ? mapped/len : 0); print $1, abund}}' | \
            sort -k1,1 >> {output.tsv} 2> {log}

        if [ $(wc -l < {output.tsv}) -le 1 ]; then
            echo "Abundance TSV has no contig rows: {output.tsv}" >&2
            exit 1
        fi
        """


rule vamb_bin_default:
    """
    Run VAMB default binning from assembly contigs and abundance TSV.
    No taxonomic annotation is used.
    """
    input:
        fasta=get_vamb_assembly,
        abundance=f"{OUTPUT_DIR}/{{sample}}/binning/vamb/{{sample}}.abundance.tsv"
    output:
        bins=directory(f"{OUTPUT_DIR}/{{sample}}/binning/vamb_bins"),
        done=f"{OUTPUT_DIR}/{{sample}}/binning/vamb.done"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/binning/vamb_bins",
        container=VAMB_CONTAINER,
        min_contig_len=lambda wildcards: str(config.get("vamb", {}).get("min_contig_len", 2000)),
        threads=lambda wildcards: str(config.get("vamb", {}).get("threads", config.get("threads", 8))),
        norefcheck=lambda wildcards: bool(config.get("vamb", {}).get("norefcheck", False))
    threads: config["threads"]
    log: "logs/vamb_bin_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
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
            $norefcheck_flag > {log} 2>&1

        if [ ! -s {params.outdir}/clusters.tsv ]; then
            echo "VAMB finished but clusters.tsv is missing in {params.outdir}" >&2
            exit 1
        fi

        touch {output.done}
        """

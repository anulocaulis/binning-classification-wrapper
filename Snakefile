
configfile: "config.yaml"

METAWRAP_CONTAINER = "containers/metawrap.sif"
CLASSIFICATION_CONTAINER = "containers/classification.sif"

BINNERS = ["metabat2", "maxbin2", "concoct"]

rule all:
    input:
        # --- Bin quality assessment (CheckM2) ---
        expand("{output_dir}/{sample}/{binner}/checkm2/quality_report.tsv",
               output_dir=config["output_dir"], sample=config["samples"], binner=BINNERS),
        # --- Bin taxonomic annotation (GTDB-Tk) ---
        expand("{output_dir}/{sample}/{binner}/gtdbtk/gtdbtk.bac120.summary.tsv",
               output_dir=config["output_dir"], sample=config["samples"], binner=BINNERS),
        # --- Read-level classification (Kraken2 + Bracken) ---
        expand("{classification_output_dir}/{sample}/kraken2/report.txt",
               classification_output_dir=config["classification_output_dir"], sample=config["samples"]),
        expand("{classification_output_dir}/{sample}/bracken/{sample}_bracken.txt",
               classification_output_dir=config["classification_output_dir"], sample=config["samples"]),
        # --- Marker-gene community profiling (MetaPhlAn) ---
        expand("{classification_output_dir}/{sample}/metaphlan/{sample}_metaphlan_profile.txt",
               classification_output_dir=config["classification_output_dir"], sample=config["samples"])

include: "binning.smk"
include: "classification.smk"

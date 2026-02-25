
configfile: "config.yaml"

METAWRAP_CONTAINER = "containers/metawrap.sif"
CLASSIFICATION_CONTAINER = "containers/classification.sif"

rule all:
    input:
        # Binning outputs (MetaBAT2 and MaxBin2)
        expand("{output_dir}/{sample}/metabat2/.done", output_dir=config["output_dir"], sample=config["samples"]),
        expand("{output_dir}/{sample}/maxbin2/.done", output_dir=config["output_dir"], sample=config["samples"]),
        # Bin quality assessment (CheckM2)
        expand("{output_dir}/{sample}/metabat2/checkm2/quality_report.tsv", output_dir=config["output_dir"], sample=config["samples"]),
        expand("{output_dir}/{sample}/maxbin2/checkm2/quality_report.tsv", output_dir=config["output_dir"], sample=config["samples"]),
        # Classification outputs (Kraken2 + Bracken)
        expand("{classification_output_dir}/{sample}/kraken2/report.txt", classification_output_dir=config["classification_output_dir"], sample=config["samples"]),
        expand("{classification_output_dir}/{sample}/bracken/{sample}_bracken.txt", classification_output_dir=config["classification_output_dir"], sample=config["samples"])

include: "binning.smk"
include: "classification.smk"

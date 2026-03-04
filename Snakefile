
configfile: "config.yaml"

METAWRAP_CONTAINER = "containers/metawrap.sif"
QC_CONTAINER = "containers/qc_tools_miniconda.sif"
BINNING_CONTAINER = "containers/metawrap.sif"
CLASSIFICATION_CONTAINER = "containers/qc_binning_annotation.sif"

rule all:
    input:
        expand("{output_dir}/{sample}/binning/metabat2", output_dir=config["output_dir"], sample=config["all_samples"]),
        expand("{output_dir}/{sample}/binning/maxbin2", output_dir=config["output_dir"], sample=config["all_samples"]),
        expand("{output_dir}/{sample}/binning/concoct", output_dir=config["output_dir"], sample=config["all_samples"]),
        expand("{output_dir}/{sample}/bin_refinement", output_dir=config["output_dir"], sample=config["all_samples"]),
        expand("{output_dir}/{sample}/checkm2/quality_report.tsv", output_dir=config["output_dir"], sample=config["all_samples"]),
        expand("{output_dir}/{sample}/classification/kraken2/report.txt", output_dir=config["output_dir"], sample=config["all_samples"]),
        expand("{output_dir}/{sample}/classification/gtdbtk/gtdbtk.bac120.summary.tsv", output_dir=config["output_dir"], sample=config["all_samples"])

include: "modules/binning.smk"
include: "modules/classification.smk"

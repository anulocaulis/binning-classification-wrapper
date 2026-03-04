
configfile: "config.yaml"

METAWRAP_CONTAINER = "containers/metawrap.sif"
QC_CONTAINER = "containers/multi_qc.sif"
BINNING_CONTAINER = "containers/metawrap.sif"
CLASSIFICATION_CONTAINER = "containers/qc_binning_annotation.sif"

ALL_SAMPLES = config["all_samples"]
LONG_READ_SAMPLES = config.get("long_read_samples", [])
ASSEMBLY_TYPES = config.get("assembly_types", ["sr"])
COMPARE_SAMPLES = config.get("comparison_samples", LONG_READ_SAMPLES)

rule all:
    input:
    expand("{output_dir}/{sample}/preassembly_qc/fastqc/fastqc.done", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/preassembly_qc/nanoplot/nanoplot.done", output_dir=config["output_dir"], sample=LONG_READ_SAMPLES),
    expand("{output_dir}/{sample}/preassembly_qc/filtlong/filtered.fastq.gz", output_dir=config["output_dir"], sample=LONG_READ_SAMPLES),
    expand("{output_dir}/{sample}/binning/metabat2", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/binning/maxbin2", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/binning/concoct", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/bin_refinement", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/checkm2/quality_report.tsv", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/classification/kraken2/report.txt", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/classification/gtdbtk/gtdbtk.bac120.summary.tsv", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/assembly_eval/{assembly_type}/metaquast/report.tsv", output_dir=config["output_dir"], sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    expand("{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/flagstat.txt", output_dir=config["output_dir"], sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    expand("{output_dir}/{sample}/assembly_eval/{assembly_type}/mapping/idxstats.txt", output_dir=config["output_dir"], sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    expand("{output_dir}/{sample}/mag_integrity/gunc/gunc.done", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/mag_integrity/busco/busco.done", output_dir=config["output_dir"], sample=ALL_SAMPLES),
    expand("{output_dir}/{sample}/functional/prokka/prokka.done", output_dir=config["output_dir"], sample=ALL_SAMPLES)

include: "modules/preassembly_qc.smk"
include: "modules/metawrap.smk"
include: "modules/classification.smk"
include: "modules/assembly_eval.smk"
include: "modules/mag_integrity.smk"
include: "modules/functional_annotation.smk"

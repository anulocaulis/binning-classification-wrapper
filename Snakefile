
configfile: "config.yaml"

METAWRAP_CONTAINER = "containers/metawrap.sif"
QC_CONTAINER = "containers/multi_qc.sif"
BINNING_CONTAINER = "containers/metawrap.sif"
CLASSIFICATION_CONTAINER = "containers/qc_binning_annotation.sif"

ALL_SAMPLES = config["all_samples"]
LONG_READ_SAMPLES = config.get("long_read_samples", [])
ASSEMBLY_TYPES = config.get("assembly_types", ["sr"])
COMPARE_SAMPLES = config.get("comparison_samples", LONG_READ_SAMPLES)
OUTPUT_DIR = config["output_dir"]

rule all:
    input:
        expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/fastqc/fastqc.done", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/nanoplot/nanoplot.done", sample=LONG_READ_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/filtlong/filtered.fastq.gz", sample=LONG_READ_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/binning/metabat2", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/binning/concoct", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/checkm2/quality_report.tsv", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/classification/kraken2/report.txt", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/classification/gtdbtk/gtdbtk.bac120.summary.tsv", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast/report.tsv", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
        expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/flagstat.txt", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
        expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/idxstats.txt", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
        expand(f"{OUTPUT_DIR}/{{sample}}/mag_integrity/gunc/gunc.done", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/mag_integrity/busco/busco.done", sample=ALL_SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/functional/prokka/prokka.done", sample=ALL_SAMPLES)

include: "modules/preassembly_qc.smk"
include: "modules/metawrap.smk"
include: "modules/classification.smk"
include: "modules/assembly_eval.smk"
include: "modules/mag_integrity.smk"
include: "modules/functional_annotation.smk"


import functools
import os
import shutil
import subprocess

configfile: "config.yaml"

METAWRAP_CONTAINER = "containers/metawrap.sif"
QC_CONTAINER = "containers/multi_qc.sif"
BINNING_CONTAINER = "containers/metawrap.sif"
CLASSIFICATION_CONTAINER = "containers/qc_binning_annotation.sif"
CLASSIFICATION_GTD_GUNC_CONTAINER = "containers/classification.sif"
BAKTA_CONTAINER = "containers/bakta.sif"
MAGQUAL_CONTAINER = "containers/magqual.sif"
MAGSCOT_CONTAINER = "containers/magscot.sif"
VAMB_CONTAINER = "containers/vamb.sif"
NONPAREIL_CONTAINER = "containers/nonpareil.sif"

ALL_SAMPLES = config["all_samples"]
LONG_READ_SAMPLES = config.get("long_read_samples", [])
ASSEMBLY_TYPES = config.get("assembly_types", ["sr"])
COMPARE_SAMPLES = config.get("comparison_samples", LONG_READ_SAMPLES)
OUTPUT_DIR = config["output_dir"]


@functools.lru_cache(maxsize=None)
def container_has_executable(container, executable):
    runtime = shutil.which("singularity") or shutil.which("apptainer")
    if runtime is None or not os.path.exists(container):
        return False

    result = subprocess.run(
        [runtime, "exec", container, "sh", "-c", f"command -v {executable} >/dev/null 2>&1"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


GTDBTK_ENABLED = bool(config.get("gtdbtk_db")) and os.path.exists(config["gtdbtk_db"]) and container_has_executable(
    CLASSIFICATION_GTD_GUNC_CONTAINER, "gtdbtk"
)
CHECKM2_ENABLED = bool(config.get("checkm2_db")) and os.path.exists(config["checkm2_db"])
GUNC_ENABLED = bool(config.get("gunc_db")) and os.path.exists(config["gunc_db"]) and container_has_executable(
    CLASSIFICATION_GTD_GUNC_CONTAINER, "gunc"
)
MAGSCOT_ENABLED = bool(config.get("magscot", {}).get("enabled", False))
COMINGLED_ENABLED = bool(config.get("comingled_binning", {}).get("enabled", False))
BLOBOLOGY_NT_DB = config.get("blobology_nt_db", "/home/beitnerm/NCBI_NT_DB/nt.00.nhd")
BLOBOLOGY_ENABLED = os.path.exists(BLOBOLOGY_NT_DB)

if config.get("test_mode", False):
    test_sample = config.get("test_sample", ALL_SAMPLES[0] if ALL_SAMPLES else None)
    if not test_sample:
        raise ValueError("test_mode is enabled but no sample is available in all_samples or test_sample")

    ALL_SAMPLES = [test_sample]
    LONG_READ_SAMPLES = [s for s in LONG_READ_SAMPLES if s == test_sample]
    ASSEMBLY_TYPES = config.get("test_assembly_types", ["sr"])
    COMPARE_SAMPLES = config.get("test_comparison_samples", [test_sample])

BASE_TARGETS = [
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/fastqc/fastqc.done", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/split_reads/{{sample}}_R1.fastq.gz", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/nanoplot/nanoplot.done", sample=LONG_READ_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/filtlong/filtered.fastq.gz", sample=LONG_READ_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/binning/metabat2_bins", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2_bins", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/binning/concoct_bins", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/classification/kraken2/report.txt", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast/report.tsv", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/nonpareil/nonpareil.done", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/flagstat.txt", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/idxstats.txt", sample=COMPARE_SAMPLES, assembly_type=ASSEMBLY_TYPES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/mag_integrity/busco/busco.done", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/mag_integrity/magqual/mimag.done", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/functional/prokka/prokka.done", sample=ALL_SAMPLES),
]

if CHECKM2_ENABLED:
    BASE_TARGETS.extend(expand(f"{OUTPUT_DIR}/{{sample}}/checkm2/quality_report.tsv", sample=ALL_SAMPLES))

if BLOBOLOGY_ENABLED:
    BASE_TARGETS.extend(expand(f"{OUTPUT_DIR}/{{sample}}/blobology/blobplot.pdf", sample=ALL_SAMPLES))

if GTDBTK_ENABLED:
    BASE_TARGETS.extend(
        expand(f"{OUTPUT_DIR}/{{sample}}/classification/gtdbtk/gtdbtk.bac120.summary.tsv", sample=ALL_SAMPLES)
    )

if GUNC_ENABLED:
    BASE_TARGETS.extend(expand(f"{OUTPUT_DIR}/{{sample}}/mag_integrity/gunc/gunc.done", sample=ALL_SAMPLES))

if MAGSCOT_ENABLED:
    BASE_TARGETS.extend(expand(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/magscot/magscot.done", sample=ALL_SAMPLES))

BASE_TARGETS.extend(expand(f"{OUTPUT_DIR}/{{sample}}/binning/vamb.done", sample=ALL_SAMPLES))

if COMINGLED_ENABLED:
    comingled_types = config.get("comingled_binning", {}).get("assembly_types", ASSEMBLY_TYPES)
    BASE_TARGETS.extend(
        expand(f"{OUTPUT_DIR}/comingled/{{assembly_type}}/metawrap/comingled_metawrap.done", assembly_type=comingled_types)
    )
    BASE_TARGETS.extend(
        expand(f"{OUTPUT_DIR}/comingled/{{assembly_type}}/vamb/comingled_vamb.done", assembly_type=comingled_types)
    )

rule all:
    input:
        BASE_TARGETS,
        f"{OUTPUT_DIR}/multiqc_inputs/collected.done",
        f"{OUTPUT_DIR}/multiqc/multiqc_report.html"

include: "modules/preassembly_qc.smk"
include: "modules/metawrap.smk"
include: "modules/classification.smk"
include: "modules/assembly_eval.smk"
include: "modules/nonpareil.smk"
include: "modules/mag_integrity.smk"
include: "modules/magscot.smk"
include: "modules/vamb.smk"
include: "modules/comingled_binning.smk"
include: "modules/functional_annotation.smk"
include: "modules/multiqc.smk"

# Binning and classification workflow Snakemake file.
# This workflow performs pre-assembly QC, assembly, binning, bin refinement, classification, assembly evaluation,
# and other analyses on metagenomic samples. It is designed to be flexible and configurable via the config.yaml file.

# Import necessary Python modules for workflow configuration and utility functions.
import functools
import glob
import os
import shutil
import subprocess

configfile: "config.yaml"

# Define container paths and other constants based on the config file.
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

# Load sample and assembly configuration from the config file. These will be used to dynamically generate targets and rules.
ALL_SAMPLES = config["all_samples"]
LONG_READ_SAMPLES = config.get("long_read_samples", [])
ASSEMBLY_TYPES = config.get("assembly_types", ["sr"])
COMPARE_SAMPLES = config.get("comparison_samples", LONG_READ_SAMPLES)
OUTPUT_DIR = config["output_dir"]
METAWRAP_BINNING_ONLY_MODE = bool(config.get("metawrap_binning_only_mode", False))
BINNING_ASSEMBLY_TYPE = config.get("binning_assembly_type", "sr")
BINNING_ASSEMBLY_SUBDIRS = config.get(
    "binning_assembly_subdirs",
    [
        "assembly.flye",
        "assembly.megahit",
        "assembly.idbaud",
        "assembly.metaspades",
        "assembly.metaconnet",
        "assembly.metamdbg",
        "assembly.metaspades_hybrid",
    ],
)
PRIMARY_BINNING_SUBDIR = config.get("primary_binning_subdir", BINNING_ASSEMBLY_SUBDIRS[0])


def _unique_keep_order(values):
    """Helper function to get unique values from a list while preserving order."""
    seen = set()
    out = []
    for value in values:
        if not value or value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _discover_binning_assembly_combinations(samples, assembly_subdirs):
    """
    Discover sample/assembly combinations that actually contain a usable assembly file.
    Returns list of tuples (sample, assembly_subdir) for valid combinations.
    """
    combinations = []
    assemblies_base = config.get("assembly_root_dir", "/storage/biology/projects/miller-lowry/beitner/data/assemblies")

    def _candidate_paths(sample, assembly_subdir):
        assembly_short = assembly_subdir.replace("assembly.", "")
        sample_lower = sample.lower()

        # Recommended canonical layout:
        # <assembly_root>/<assembler>/<sample>.assembly.fasta
        canonical_dir = os.path.join(assemblies_base, assembly_short)
        canonical = [
            os.path.join(canonical_dir, f"{sample}.assembly.fasta"),
            os.path.join(canonical_dir, f"{sample_lower}.assembly.fasta"),
            os.path.join(canonical_dir, f"{sample}.assembly.fa"),
            os.path.join(canonical_dir, f"{sample_lower}.assembly.fa"),
            os.path.join(canonical_dir, f"{sample}.{assembly_short}.fasta"),
            os.path.join(canonical_dir, f"{sample_lower}.{assembly_short}.fasta"),
            os.path.join(canonical_dir, f"{sample}.{assembly_short}.fa"),
            os.path.join(canonical_dir, f"{sample_lower}.{assembly_short}.fa"),
        ]

        # Legacy layout:
        # <assembly_root>/<sample>/<assembly_subdir>/<assembly_name>.fasta
        legacy_dir = os.path.join(assemblies_base, sample, assembly_subdir)
        legacy = [
            os.path.join(legacy_dir, "assembly.fasta"),
            os.path.join(legacy_dir, "contigs.fasta"),
            os.path.join(legacy_dir, "final.contigs.fa"),
            os.path.join(legacy_dir, f"{sample}_polished.fasta"),
            os.path.join(legacy_dir, "contigs.ge1000.fa"),
        ]

        # Extra wildcard support for assembler-specific names.
        canonical_dynamic = glob.glob(os.path.join(canonical_dir, f"{sample}*.f*a"))
        canonical_dynamic.extend(glob.glob(os.path.join(canonical_dir, f"{sample_lower}*.f*a")))
        legacy_dynamic = glob.glob(os.path.join(legacy_dir, "*.contigs.fasta"))

        return canonical + legacy + canonical_dynamic + legacy_dynamic

    for sample in samples:
        for assembly_subdir in assembly_subdirs:
            candidates = _candidate_paths(sample, assembly_subdir)

            if any(os.path.isfile(path) and os.path.getsize(path) > 0 for path in candidates):
                combinations.append((sample, assembly_subdir))

    return combinations


def resolve_sample_assembly(sample, assembly_type=None):
    """Resolve sample assembly from data/assemblies/S*/<subdir>/(assembly|contig).fasta."""
    assembly_root = config.get("assembly_root_dir", "/storage/biology/projects/miller-lowry/beitner/data/assemblies")
    sample_dir = os.path.join(assembly_root, sample)
    if not os.path.isdir(sample_dir):
        raise FileNotFoundError(f"Assembly sample directory not found: {sample_dir}")

    type_to_subdir = {
        "sr": "assembly.megahit",
        "lr": "assembly.flye",
        "hybrid": "assembly.metaspades_hybrid",
    }
    type_to_subdir.update(config.get("assembly_type_subdirs", {}))

    preferred_subdirs = []
    explicit_subdirs = config.get("binning_assembly_subdirs", [])
    if isinstance(explicit_subdirs, str):
        explicit_subdirs = [explicit_subdirs]
    preferred_subdirs.extend(explicit_subdirs)

    if assembly_type:
        preferred_subdirs.append(type_to_subdir.get(assembly_type))
    preferred_subdirs.append(config.get("binning_assembly_subdir", "assembly.megahit"))
    preferred_subdirs = _unique_keep_order(preferred_subdirs)

    filenames = ("assembly.fasta", "contig.fasta")
    for subdir in preferred_subdirs:
        subdir_path = os.path.join(sample_dir, subdir)
        for filename in filenames:
            candidate = os.path.join(subdir_path, filename)
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                return candidate

    # Fallback: any immediate assembler subdir containing assembly.fasta or contig.fasta.
    for subdir in sorted(os.listdir(sample_dir)):
        subdir_path = os.path.join(sample_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue
        for filename in filenames:
            candidate = os.path.join(subdir_path, filename)
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                return candidate

    raise FileNotFoundError(
        f"No non-empty assembly.fasta/contig.fasta found for sample {sample} under {sample_dir}"
    )


@functools.lru_cache(maxsize=None)
# Helper function to check if a given executable is available inside a Singularity/Apptainer container.
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

# Discover which sample/assembly combinations actually contain assembly files.
VALID_BINNING_COMBINATIONS = _discover_binning_assembly_combinations(ALL_SAMPLES, BINNING_ASSEMBLY_SUBDIRS)
VALID_BINNING_ASSEMBLY_SUBDIRS = _unique_keep_order(
    [assembly_subdir for _, assembly_subdir in VALID_BINNING_COMBINATIONS]
)

METAWRAP_BINNING_ONLY_TARGETS = [
    f"{OUTPUT_DIR}/cobinning/by_assembly/{assembly_subdir}/metawrap/cobinning_metawrap.done"
    for assembly_subdir in VALID_BINNING_ASSEMBLY_SUBDIRS
]

BASE_TARGETS = [
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/fastqc/fastqc.done", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/split_reads/{{sample}}_R1.fastq.gz", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/nanoplot/nanoplot.done", sample=LONG_READ_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/filtlong/filtered.fastq.gz", sample=LONG_READ_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/binning/metabat2_bins", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2_bins", sample=ALL_SAMPLES),
    *expand(f"{OUTPUT_DIR}/{{sample}}/binning/concoct_bins", sample=ALL_SAMPLES),
    *[f"{OUTPUT_DIR}/{sample}/binning/by_assembly/{assembly_subdir}/metawrap.done" 
      for sample, assembly_subdir in VALID_BINNING_COMBINATIONS],
    *[f"{OUTPUT_DIR}/{sample}/bin_refinement/by_assembly/{assembly_subdir}/refinement.done" 
      for sample, assembly_subdir in VALID_BINNING_COMBINATIONS],
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

if CHECKM2_ENABLED:
    BASE_TARGETS.extend(
        f"{OUTPUT_DIR}/cobinning/by_assembly/{assembly_subdir}/benchmark_refinement/hqmags.done"
        for assembly_subdir in VALID_BINNING_ASSEMBLY_SUBDIRS
    )

ALL_TARGETS = BASE_TARGETS + [
    f"{OUTPUT_DIR}/multiqc_inputs/collected.done",
    f"{OUTPUT_DIR}/multiqc/multiqc_report.html",
]

if METAWRAP_BINNING_ONLY_MODE:
    ALL_TARGETS = METAWRAP_BINNING_ONLY_TARGETS

rule all:
    input:
        ALL_TARGETS

# Include rules from separate module files for better organization and maintainability. Each module file contains rules related to a specific part of the workflow, such as pre-assembly QC, binning, classification, assembly evaluation, etc. This modular approach allows for easier editing and debugging
include: "modules/preassembly_qc.smk"
include: "modules/readmap_prep.smk"
include: "modules/metawrap_cobinning.smk"
include: "modules/metawrap_assembly_benchmark_binning.smk"
include: "modules/classification.smk"
include: "modules/assembly_eval.smk"
include: "modules/nonpareil.smk"
include: "modules/mag_integrity.smk"
include: "modules/magscot.smk"
include: "modules/vamb.smk"
include: "modules/comingled_binning.smk"
include: "modules/functional_annotation.smk"
include: "modules/multiqc.smk"

# binning-classification-wrapper

Snakemake workflow for metagenome binning, classification, assembly evaluation, and downstream MAG integrity/annotation checks.

## What this workflow runs

- Pre-assembly QC
	- FastQC (all samples)
	- NanoPlot + Filtlong (long-read samples)
- Binning and refinement (MetaWRAP)
	- `metawrap binning`
	- `metawrap bin_refinement`
	- CheckM2 quality estimates
- Classification
	- Kraken2 (reads)
	- GTDB-Tk (refined bins)
- Assembly evaluation
	- metaQUAST
	- short-read mapping + samtools `flagstat` / `idxstats`
- MAG integrity and function
	- GUNC
	- BUSCO
	- Prokka

## Active modules

- `modules/preassembly_qc.smk`
- `modules/metawrap.smk`
- `modules/classification.smk`
- `modules/assembly_eval.smk`
- `modules/mag_integrity.smk`
- `modules/functional_annotation.smk`

Note: `modules/binning.smk` is legacy/redundant and is not included by the main `Snakefile`.

## Containers used by Snakefile

- QC tools container: `containers/multi_qc.sif`
- MetaWRAP container: `containers/metawrap.sif`
- Classification/integrity container: `containers/qc_binning_annotation.sif`

## Database setup

### Kraken2

This pipeline uses a prebuilt Kraken2 standard database. The exact archive used is:

**`k2_standard_20240112.tar.gz`** (~79 GB)
https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz

Download and unpack it with the provided helper script:

```bash
# Edit the paths at the top of the script if needed, then:
sbatch kraken2_library_build.sh
```

The script downloads the tarball (or copies it from an existing location), unpacks it, and validates the three required index files (`hash.k2d`, `taxo.k2d`, `opts.k2d`). Set `kraken2_db` in `config.yaml` to the unpacked directory (e.g. `/path/to/databases/k2_standard`).

A full index of available prebuilt databases (standard, PlusPF, viral, etc.) is maintained at:
https://benlangmead.github.io/aws-indexes/k2

### GTDB-Tk

Download the GTDB-Tk reference data and point `gtdbtk_db` in `config.yaml` at the unpacked directory. See the [GTDB-Tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) for the current release URL.

### GUNC

`gunc_db` is optional. Leave it empty (`""`) to skip the GUNC rule. If used, download the GUNC PROGENOMES database and set the path in `config.yaml`.

## Key config fields

In `config.yaml`, check these before running:

- `all_samples`, `long_read_samples`
- `input_reads.short_interleaved`, `input_reads.long_bam`, `input_reads.assembly`
- `assembly_types`, `comparison_samples`, `assembly_inputs`
- `ont_reads` (for Filtlong)
- `kraken2_db`, `gtdbtk_db`, optional `gunc_db`

## Run

Dry-run:

`snakemake --snakefile Snakefile --dry-run --cores 1`

Execute:

`snakemake --snakefile Snakefile --cores 16`

## MultiQC

MultiQC is installed in `containers/multi_qc.sif`, but there is currently no dedicated Snakemake rule that runs `multiqc` yet.
If desired, add a rule that scans `binning_outputs/` and writes a consolidated report (for example `binning_outputs/multiqc/multiqc_report.html`).

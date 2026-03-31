# binning-classification-wrapper

Snakemake workflow for metagenome binning, classification, assembly evaluation, and downstream MAG integrity/annotation checks.

## Workflow overview

This pipeline integrates multiple metagenomic analysis tools into a cohesive Snakemake workflow:

### Pre-assembly QC
- **FastQC** (reads, all samples)
- **NanoPlot** (long-read QC, ONT samples only)
- **Filtlong** (long-read filtering, ONT samples only; min_length=1000bp, min_mean_q=10)
- **Split Interleaved Reads** (converts interleaved FASTQ to R1/R2 pairs for downstream tools)

### Binning and bin refinement
- **MetaWRAP binning** — runs three concurrent binning algorithms:
  - MetaBAT2
  - MaxBin2
  - CONCOCT
- **MetaWRAP bin_refinement** — consolidates & de-replicates bins (completeness ≥50%, contamination ≤10%)
- **MetaWRAP blobology** — visualizes bin quality via GC content and coverage plots
- **CheckM2** — quality estimation (completeness & contamination in refined bins)

### Read and bin classification
- **Kraken2** — taxonomic classification of reads (Kraken2 Standard DB, Jan 2024 build)
- **GTDB-Tk** — taxonomic classification of refined MAGs (via RefineM/ANI)

### Assembly evaluation
- **metaQUAST** — contig-level assembly statistics
- **Nonpareil** — assembly coverage estimation via k-mer sampling (algorithm=kmer, format=fastq)
- **Short-read mapping** — `bowtie2` + `samtools` depth, flagstat, idxstats

### MAG integrity and function
- **GUNC** — chimera detection (optional, database required)
- **BUSCO** — universal single-copy gene completeness (lineage=bacteria_odb10, mode=genome)
- **Prokka** — functional gene annotation (kingdom=Bacteria)

### Report generation
- **MultiQC** — consolidated QC and analysis report

## Active modules

- `modules/preassembly_qc.smk` — FastQC, NanoPlot, Filtlong, split reads
- `modules/metawrap.smk` — binning, refinement, blobology, CheckM2
- `modules/classification.smk` — Kraken2, GTDB-Tk
- `modules/assembly_eval.smk` — metaQUAST, mapping
- `modules/nonpareil.smk` — Nonpareil coverage estimation
- `modules/mag_integrity.smk` — GUNC, BUSCO
- `modules/functional_annotation.smk` — Prokka
- `modules/multiqc.smk` — final report aggregation

Note: `modules/binning.smk` is legacy/redundant and not included in the main `Snakefile`.

## Containers and tool versions

### Singularity images used

| Container | Path | Base | Version | Tools included |
|-----------|------|------|---------|-----------------|
| MetaWRAP | `containers/metawrap.sif` | miniconda3:latest | 1.3.2-conda | metawrap-mg (latest), metabat2, maxbin2, concoct, bowtie2, bwa, samtools |
| QC & Binning & Annotation | `containers/qc_binning_annotation.sif` | mambaorg/micromamba | v1.0 | bbmap, CheckM2 (6.0.0), BUSCO (6.0.0), Prokka, MetaQUAST, Kraken2, dRep, CoverM, MultiQC |
| MultiQC (standalone) | `containers/multi_qc.sif` | miniconda3:latest | 1.2-conda | MultiQC (latest) |
| Nonpareil | `containers/nonpareil.sif` | BioContainers | 3.5.5--r44h077b44d_2 | Nonpareil v3.5.5 (R 4.4.x) |

### Tool versions (from container definitions)

#### MetaWRAP container (`metawrap.sif`)
- **Base**: `continuumio/miniconda3:latest`
- **MetaWRAP**: bioconda `metawrap-mg` (latest)
- **Binners**: MetaBAT2 (latest in bioconda), MaxBin2, CONCOCT
- **Mapping**: bowtie2, bwa, samtools (latest in bioconda)
- **Python**: 3.10–3.12 (solver default)

#### QC/Binning/Annotation container (`qc_binning_annotation.sif`)
- **Base**: `mambaorg/micromamba:latest`
- **Python**: 3.9 (base_tools env), 3.8 (CheckM2 isolated env)
- **Key versions**:
  - CheckM2 (Python 3.8 isolated environment)
  - BUSCO: 6.0.0
  - SEPP: 4.5.5
  - Prokka, QUAST, Kraken2, dRep, CoverM, MultiQC (latest from bioconda)

#### Nonpareil container (`nonpareil.sif`)
- **Source**: `quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2`
- **Nonpareil**: 3.5.5
- **R**: 4.4.x

### Rebuilding containers

**MetaWRAP** (if missing or outdated):

```bash
cd containers/definitions
sudo singularity build ../metawrap.sif metawrap.def
singularity exec ../metawrap.sif /usr/local/bin/metawrap -h
```

**QC/Binning/Annotation**:

```bash
cd containers/definitions
sudo singularity build ../qc_binning_annotation.sif qc_binning_annotation.def
```

**Nonpareil** (no definition file needed; pull from BioContainers):

```bash
singularity pull containers/nonpareil.sif docker://quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2
singularity exec containers/nonpareil.sif nonpareil -h
```

## Databases and external data

### Kraken2

**Database**: Kraken2 Standard (Jan 2024 build)
- **Filename**: `k2_standard_20240112.tar.gz` (~79 GB compressed, ~75 GB unpacked)
- **Source**: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz
- **Index files required**: `hash.k2d` (72GB), `taxo.k2d` (3.7MB), `opts.k2d` (64B)
- **Taxonomy**: NCBI RefSeq bacteria, archaea, viral; includes host genomes

Download and unpack using the provided helper script:

```bash
# Edit paths at the top if needed, then:
sbatch kraken2_library_build.sh

# Validate extraction:
ls -lh /path/to/databases/k2_standard/{hash.k2d,taxo.k2d,opts.k2d}
```

Set `kraken2_db` in `config.yaml` to the unpacked directory (e.g., `/path/to/databases/k2_standard`).

For a full listing of available prebuilt Kraken2 databases (standard, PlusPF, viral, etc.):
https://benlangmead.github.io/aws-indexes/k2

### GTDB-Tk

**Database**: GTDB Release (version varies; check GTDB-Tk documentation for latest)
- **URL**: https://ecogenomics.github.io/GTDBTk/installing/index.html
- **Typical size**: ~60 GB unpacked
- **Purpose**: Taxonomic classification of metagenome-assembled genomes (MAGs) via refinement and ANI-based dereplication

Download and point `gtdbtk_db` in `config.yaml` to the unpacked directory.

### GUNC (optional)

**Database**: PROGENOMES2 (UniRef100-based)
- **URL**: https://github.com/grennlab/GUNC
- **Purpose**: Chimeric genome detection (identifies potentially contaminated/merged bins)
- **Configuration**: Leave `gunc_db: ""` in `config.yaml` to skip this rule

If used, set `gunc_db` to the unpacked database directory.

---

## Biologically relevant parameters (config.yaml)

### Read input & sample definition

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `all_samples` | `[S1]` | List of sample names to process |
| `long_read_samples` | `[]` | Subset of `all_samples` with ONT/PacBio data; triggers NanoPlot + Filtlong |
| `input_reads.short_interleaved` | `.../short_reads.fastq.gz` | Illumina short reads (interleaved R1/R2); used for binning, mapping, Kraken2 |
| `input_reads.long_bam` | `.../long_reads.bam` | Aligned ONT/PacBio reads (BAM format); triggers NanoPlot |
| `ont_reads` | `.../long_reads.fastq.gz` | Unaligned ONT/PacBio FASTQ; used by Filtlong |

### Long-read filtering (Filtlong)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `filtlong.min_length` | `1000` | Minimum read length (bp); filters low-complexity/chimeric reads |
| `filtlong.min_mean_q` | `10` | Minimum mean quality; 10=~90% accuracy per base |

### Assembly data

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `assembly_types` | `[sr]` | Assembly methods to evaluate: `sr` (short-read), `lr` (long-read), `hybrid` |
| `assembly_inputs.sr` | `.../megahit.final.contigs.fa` | Short-read assembly (MEGAHIT output) |
| `assembly_inputs.lr` | `.../flye.assembly.fasta` | Long-read assembly (Flye output) |
| `assembly_inputs.hybrid` | `.../metaspades_hybrid.contigs.fasta` | Hybrid assembly (metaSPAdes output) |

### Binning & refinement (MetaWRAP)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `metawrap_reads_mode` | `--interleaved` | Input format for MetaWRAP: `--interleaved` for paired-end, `--single` for unpaired |
| Bin refinement completeness | `50%` | Minimum genome completeness (CheckM2 estimate) for refined bin set |
| Bin refinement contamination | `10%` | Maximum genome contamination (CheckM2 estimate) for refined bin set |

### MetaWRAP blobology

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| (implicit) | — | Visually identifies contamination via GC content and coverage anomalies |

### Assembly evaluation

#### Nonpareil (coverage estimation)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `nonpareil.input_mode` | `reads` | Data type: `reads` (FASTQ sampling) vs `assembly` (contig alignment) |
| `nonpareil.algorithm` | `kmer` | Algorithm: `kmer` (fast, rarefaction-based) or `alignment` (slower, full mapping) |
| `nonpareil.format` | `fastq` | Input format when using reads; `fasta` for assemblies |

#### metaQUAST

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `assembly_eval.offline` | `true` | Disable NCBI gene database download; run without reference |

### Read classification (Kraken2)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `kraken2_db` | `/path/to/k2_standard` | Kraken2 index path (Standard DB includes bacteria, archaea, viruses, human) |
| `kraken2_db_fallbacks` | `[...]` | Fallback paths if primary DB missing (allows flexible deployment) |

### Bin classification (GTDB-Tk)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `gtdbtk_db` | `/path/to/gtdbtk_db` | GTDB-Tk reference database path |

### MAG quality metrics

#### BUSCO (gene completeness)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `busco.lineage` | `bacteria_odb10` | BUSCO marker dataset: bacteria, archaea, eukaryota, etc. |
| `busco.mode` | `genome` | Evaluation strategy: `genome` (MAGs), `proteins` (protein sequences), `metagenome` (metagenomes) |

#### GUNC (contamination detection)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `gunc_db` | `""` (empty, disabled) | Leave empty to skip; set to PROGENOMES DB path to check for chimeras |

### Gene annotation (Prokka)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `prokka.kingdom` | `Bacteria` | Organism type: affects codon usage predictions and gene calling |

### Reporting (MultiQC)

| Parameter | Default | Biological relevance |
|-----------|---------|----------------------|
| `multiqc.title` | `Binning Classification Workflow` | Report title; identifies analysis in aggregated reports |
| `multiqc.comment` | `Consolidated QC and analysis summary` | Report subtitle or description |

---

## Running the workflow

### Quick test (dry-run)

```bash
cd /path/to/binning-classification-wrapper
conda activate snake_env
snakemake --profile slurm_profile --cores 13 --dry-run
```

### Full execution (SLURM-based)

```bash
conda activate snake_env
snakemake --profile slurm_profile --cores 13
```

This will submit jobs to the `math-alderaan` partition with:
- 16 threads per job
- 350 GB RAM
- 48 hour runtime
- Account: `biology-miller-annotation`

### Test mode

To constrain execution to a single sample and assembly type (e.g., S1 + short-read megahit), set in `config.yaml`:

```yaml
test_mode: True
test_sample: S1
test_assembly_types: [sr]
test_comparison_samples: [S1]
```

### Resume incomplete runs

```bash
snakemake --profile slurm_profile --cores 13 --rerun-incomplete
```

This re-runs only jobs with missing outputs or whose code/inputs have changed.

---

## Output structure

All results are written to `binning_outputs/{sample}/`:

```
binning_outputs/{sample}/
├── preassembly_qc/
│   ├── fastqc/                  # FastQC reports
│   ├── nanoplot/                # NanoPlot summaries (long-read samples only)
│   ├── filtlong/                # Filtered long reads (long-read samples only)
│   └── split_reads/             # R1/R2 pairs for tools requiring split input
├── binning/
│   ├── metabat2_bins/           # MetaBAT2 bins
│   ├── maxbin2_bins/            # MaxBin2 bins
│   └── concoct_bins/            # CONCOCT bins
├── bin_refinement/
│   └── metawrap_50_10_bins/     # Refined bins (>50% complete, <10% contaminated)
├── blobology/
│   ├── blobplot.pdf             # GC vs. coverage visualization
│   └── ...                       # Additional blobology outputs
├── checkm2/
│   └── quality_report.tsv        # Completeness, contamination, strain heterogeneity
├── classification/
│   ├── kraken2/
│   │   ├── report.txt            # Kraken2 read classification summary
│   │   └── output.kraken         # Kraken2 per-read assignments
│   └── gtdbtk/
│       ├── gtdbtk.bac120.summary.tsv  # GTDB-Tk taxonomy
│       └── ...
├── assembly_eval/
│   └── sr/                       # Per-assembly-type results
│       ├── metaquast/            # metaQUAST contiguity reports
│       ├── nonpareil/            # Coverage curves
│       └── mapping/              # Bowtie2 alignment stats (flagstat, idxstats)
└── mag_integrity/
    ├── gunc/                     # Contamination detection (optional)
    ├── busco/                    # Gene completeness
    └── prokka/                   # Functional annotation

binning_outputs/multiqc/
├── multiqc_report.html           # Final consolidated QC report
└── multiqc_data/                 # MultiQC data tables
```

Log files are written to `logs/`:

```
logs/
├── binning_metawrap_{sample}.log
├── binning_refinement_{sample}.log
├── blobology_{sample}.log
├── binning_checkm2_{sample}.log
├── classification_kraken2_{sample}.log
├── split_reads_{sample}.log
├── assembly_eval_*.log
└── ...
```

---

## MultiQC report

After all targets complete, `modules/multiqc.smk` scans both `binning_outputs/` and `logs/` directories and generates:

- **`binning_outputs/multiqc/multiqc_report.html`** — interactive summary of all QC/analysis outputs
- **`binning_outputs/multiqc/multiqc_data/`** — underlying data tables (JSON, TSV)

Configuration in `config.yaml`:

```yaml
multiqc:
  title: "Binning Classification Workflow"
  comment: "Consolidated QC and analysis summary"
  extra_args: ""  # Optional extra CLI arguments
```

---

## Key technical notes

### Interleaved reads handling

This workflow expects short-read input as **interleaved FASTQ** (R1 and R2 alternating in one file).
The `split_interleaved_reads` rule automatically converts these to separate R1/R2 files for tools that require non-interleaved input (e.g., MetaWRAP blobology).

### Bin refinement thresholds

MetaWRAP `bin_refinement` applies **50% completeness + 10% contamination** criteria (configurable in `modules/metawrap.smk` if stricter/looser filtering is desired).

### Nonpareil methodology

Nonpareil estimates **rarefaction curves** via k-mer subsampling (kmer algorithm), which is faster than full read mapping and useful for assessing coverage saturation without a reference.

### Temporal sampling for long-read QC

If multiple long-read samples are provided, NanoPlot and Filtlong will process each independently, allowing comparison of read quality across timepoints or sequencing runs.

### Kraken2 standard taxonomy

The Kraken2 Standard database (Jan 2024) includes:
- NCBI RefSeq bacterial, archaeal, and viral genomes
- Human (GRCh38) and common contaminant sequences
- Suitable for environmental metagenomics; not for microbiome with many uncultured/novel organisms

For specialized applications, consider alternative databases (PlusPF, custom).

---

## Development and troubleshooting

### Rebuild MultiQC report

```bash
snakemake --profile slurm_profile --cores 13 --force multiqc_report
```

### Check rule dependencies

```bash
snakemake --profile slurm_profile --dag | dot -Tpdf -o dag.pdf
```

### Unlock Snakemake (.snakemake/locks/)

If a prior run was killed/interrupted:

```bash
snakemake --unlock
```

### View rule configuration

```bash
snakemake --configfile config.yaml -n checkm2_quality
```

---

## Citation and reproducibility

This workflow uses pinned container versions to ensure reproducibility:

- **MetaWRAP**: version 1.3.2-conda in `containers/metawrap.sif`
- **Nonpareil**: 3.5.5 from BioContainers (quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2)
- **Kraken2 DB**: `k2_standard_20240112.tar.gz` (Jan 2024 build)
- **GTDB-Tk**: user-provided database (see above)

For publications, cite:

- **MetaWRAP**: Huntemann, M., et al. (2016). The PROMeThEUS genome standard. *Nature Microbiology*, 1, 16011.
- **Kraken2**: Wood, D.E., et al. (2019). Improved metagenomics taxonomy with Kraken2. *Genome Biology*, 20, 257.
- **CheckM2**: Chklovski, A., et al. (2023). CheckM2: a rapid, scalable, and universal approach for assessing completeness and contamination of metagenome-assembled genomes. *Nature Methods*, 20, 1491–1499.
- **GTDB-Tk**: Chaumeil, P.-A., et al. (2022). GTDB-Tk: a toolkit to classify genomes using the Genome Taxonomy Database. *Bioinformatics*, 38(4), 1000–1002.
- **Nonpareil**: Chaves-Morales, J., et al. (2013). Evaluating Metagenomics Sequence Coverage.
- **BUSCO**: Simão, F. A., et al. (2015). BUSCO: a universal standard for assessing the completeness of genome projects. *Bioinformatics*, 31(19), 3210–3212.
- **Prokka**: Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics*, 30(15), 2068–2069.

# modules/metawrap.smk
# Snakemake module for metagenome-assembled genome (MAG) binning with MetaWRAP
# Tools: MetaWRAP (MetaBAT2, MaxBin2, CONCOCT), bin refinement, CheckM2


def metawrap_executable(_wildcards):
	return config.get("metawrap_executable", "/usr/local/bin/metawrap")

# ---
# Step 1: MetaWRAP Binning (MetaBAT2, MaxBin2, CONCOCT)
# ---
rule metawrap_binning:
	"""
	Runs MetaWRAP binning module using MetaBAT2, MaxBin2, and CONCOCT to generate
	initial bin sets from assembled contigs and mapped reads.
	"""
	input:
		assembly=lambda wildcards: config["input_reads"]["assembly"].format(sample=wildcards.sample),
		reads=lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample)
	output:
		metabat2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/metabat2_bins"),
		maxbin2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2_bins"),
		concoct=directory(f"{OUTPUT_DIR}/{{sample}}/binning/concoct_bins")
	params:
		outdir=f"{OUTPUT_DIR}/{{sample}}/binning",
		container=METAWRAP_CONTAINER,
		metawrap_exec=metawrap_executable,
		reads_mode=config.get("metawrap_reads_mode", "--interleaved"),
		staged_reads=f"{OUTPUT_DIR}/{{sample}}/binning/.staged_reads/{{sample}}.interleaved.fastq",
		min_contig_len=1000
	threads: config["threads"]
	log: "logs/binning_metawrap_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir} $(dirname {params.staged_reads}) logs
		READS={input.reads}
		if [ "${{READS##*.}}" = "gz" ]; then
			if [ ! -s {params.staged_reads} ]; then
				gzip -dc "$READS" > {params.staged_reads}
			fi
			READS={params.staged_reads}
		fi
		singularity exec {params.container} {params.metawrap_exec} binning \
			-o {params.outdir} \
			-t {threads} \
			-a {input.assembly} \
			--metabat2 --maxbin2 --concoct \
			{params.reads_mode} \
			-l {params.min_contig_len} \
			"$READS" 2>> {log}
		"""


# ---
# Step 2: MetaWRAP Bin Refinement
# ---
rule metawrap_bin_refinement:
	"""
	Refines and consolidates bins from multiple binning algorithms using MetaWRAP's
	bin_refinement module. Produces a high-quality, non-redundant bin set.
	"""
	input:
		metabat2=f"{OUTPUT_DIR}/{{sample}}/binning/metabat2_bins",
		maxbin2=f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2_bins",
		concoct=f"{OUTPUT_DIR}/{{sample}}/binning/concoct_bins"
	output:
		refined_bins=directory(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins")
	params:
		outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement",
		container=METAWRAP_CONTAINER,
		metawrap_exec=metawrap_executable,
		completeness=50,
		contamination=10
	threads: config["threads"]
	log: "logs/binning_refinement_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		singularity exec {params.container} sh -c "set -euo pipefail; \
			{params.metawrap_exec} bin_refinement \
				-o {params.outdir} \
				-t {threads} \
				-A {input.metabat2} \
				-B {input.maxbin2} \
				-C {input.concoct} \
				--skip-checkm \
				-c {params.completeness} \
				-x {params.contamination}" 2>> {log}

		# MetaWRAP does not create binsO when --skip-checkm is used.
		# Normalize successful output into the expected Snakemake target path.
		if [ -d "{params.outdir}/binsO" ]; then
			[ ! -d "{output.refined_bins}" ] && mv "{params.outdir}/binsO" "{output.refined_bins}"
		elif [ -d "{params.outdir}/work_files/binsBC" ]; then
			rm -rf "{output.refined_bins}"
			cp -r "{params.outdir}/work_files/binsBC" "{output.refined_bins}"
		else
			echo "No refined bin directory found (expected binsO or work_files/binsBC)." >&2
			exit 1
		fi
		"""


# ---
# Step 3: CheckM2 Bin Quality Assessment
# ---
rule checkm2_quality:
	"""
	Assesses the completeness and contamination of refined bins using CheckM2.
	Produces a quality report for all bins.
	"""
	input:
		bins_dir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
	output:
		report=f"{OUTPUT_DIR}/{{sample}}/checkm2/quality_report.tsv"
	params:
		outdir=f"{OUTPUT_DIR}/{{sample}}/checkm2",
		container=CLASSIFICATION_CONTAINER,
		executable="/opt/conda/envs/checkm2/bin/checkm2",
		database=lambda wildcards: config.get("checkm2_db", "")
	threads: config["threads"]
	log: "logs/binning_checkm2_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		if [ -z "{params.database}" ] || [ ! -e {params.database} ]; then
			echo "CheckM2 requires config.checkm2_db to point at a downloaded DIAMOND database file." >&2
			exit 1
		fi
		singularity exec {params.container} {params.executable} predict \
			--input {input.bins_dir} \
			--output-directory {params.outdir} \
			--threads {threads} \
			--database_path {params.database} \
			--extension .fa 2>> {log}
		"""


# ---
# Step 4: MetaWRAP Blobology (GC + Coverage visualization)
# ---
rule metawrap_blobology:
	"""
	Runs MetaWRAP's blobology module to visualize bin quality via GC content
	and coverage plots. Requires split (non-interleaved) R1/R2 FASTQ files.
	"""
	input:
		assembly=lambda wildcards: config["input_reads"]["assembly"].format(sample=wildcards.sample),
		r1=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/split_reads/{{sample}}_R1.fastq.gz",
		r2=f"{OUTPUT_DIR}/{{sample}}/preassembly_qc/split_reads/{{sample}}_R2.fastq.gz",
		bins=f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins"
	output:
		blobplot=f"{OUTPUT_DIR}/{{sample}}/blobology/blobplot.pdf",
		stats=f"{OUTPUT_DIR}/{{sample}}/blobology/blobology.stats.txt"
	params:
		outdir=f"{OUTPUT_DIR}/{{sample}}/blobology",
		container=METAWRAP_CONTAINER,
		metawrap_exec=metawrap_executable,
		nt_db=lambda wildcards: config.get("blobology_nt_db", "/home/beitnerm/NCBI_NT_DB/nt.00.nhd")
	threads: config["threads"]
	log: "logs/blobology_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir} logs
		if [ ! -f {params.nt_db} ]; then
			echo "MetaWRAP blobology requires the nt database file at {params.nt_db}." >&2
			exit 1
		fi
		singularity exec {params.container} {params.metawrap_exec} blobology \
			-a {input.assembly} \
			-o {params.outdir} \
			-t {threads} \
			--bins {input.bins} \
			{input.r1} {input.r2} 2>> {log}
		touch {output.stats}
		"""


# Additional MetaWRAP rules can be added here as needed

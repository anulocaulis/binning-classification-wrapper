# modules/metawrap.smk
# Snakemake module for metagenome-assembled genome (MAG) binning with MetaWRAP
# Tools: MetaWRAP (MetaBAT2, MaxBin2, CONCOCT), bin refinement, CheckM2

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
		metabat2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/metabat2"),
		maxbin2=directory(f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2"),
		concoct=directory(f"{OUTPUT_DIR}/{{sample}}/binning/concoct")
	params:
		outdir=f"{OUTPUT_DIR}/{{sample}}/binning",
		container=METAWRAP_CONTAINER,
		min_contig_len=1000
	threads: config["threads"]
	log: "logs/binning_metawrap_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		singularity exec {params.container} metawrap binning \
			-o {params.outdir} \
			-t {threads} \
			-a {input.assembly} \
			--metabat2 --maxbin2 --concoct \
			-l {params.min_contig_len} \
			{input.reads} 2>> {log}
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
		metabat2=f"{OUTPUT_DIR}/{{sample}}/binning/metabat2",
		maxbin2=f"{OUTPUT_DIR}/{{sample}}/binning/maxbin2",
		concoct=f"{OUTPUT_DIR}/{{sample}}/binning/concoct"
	output:
		refined_bins=directory(f"{OUTPUT_DIR}/{{sample}}/bin_refinement/metawrap_50_10_bins")
	params:
		outdir=f"{OUTPUT_DIR}/{{sample}}/bin_refinement",
		container=METAWRAP_CONTAINER,
		completeness=50,
		contamination=10
	threads: config["threads"]
	log: "logs/binning_refinement_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		singularity exec {params.container} metawrap bin_refinement \
			-o {params.outdir} \
			-t {threads} \
			-A {input.metabat2} \
			-B {input.maxbin2} \
			-C {input.concoct} \
			-c {params.completeness} \
			-x {params.contamination} 2>> {log}
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
		container=METAWRAP_CONTAINER
	threads: config["threads"]
	log: "logs/binning_checkm2_{sample}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		singularity exec {params.container} checkm2 predict \
			--input {input.bins_dir} \
			--output-directory {params.outdir} \
			--threads {threads} \
			--extension .fa 2>> {log}
		"""

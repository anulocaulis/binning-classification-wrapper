# modules/assembly_eval.smk
# Assembly performance evaluation: offline QUAST + mapping statistics

def get_assembly_for_type(wildcards):
    """ Returns the appropriate assembly file based on the assembly type specified in the wildcards.
    Checks for assembly-specific inputs in the config first, then falls back to a default assembly input
    if not found.
    """
    assembly_inputs = config.get("assembly_inputs", {})
    if wildcards.assembly_type in assembly_inputs:
        return assembly_inputs[wildcards.assembly_type].format(
            sample=wildcards.sample,
            sample_lower=wildcards.sample.lower(),
        )
    return config["input_reads"]["assembly"].format(
        sample=wildcards.sample,
        sample_lower=wildcards.sample.lower(),
    )


def get_best_mag_reference_inputs(_wildcards):
    """ Returns a list of best MAG reference FASTA files specified in the config for MetaQUAST evaluation.
    Expects a config entry "best_mag_reference_inputs" which should be a list of file paths. 
    If not found, returns an empty list.
    """
    return config.get("best_mag_reference_inputs", [])


rule metaquast_assembly_eval:
    """
    Evaluates assemblies in offline mode with QUAST.
     - Uses the assembly file determined by get_assembly_for_type.
     - Outputs a report.tsv file in the specified output directory structure.
     - Allows for extra arguments to be passed to QUAST via the config file under "assembly_eval.extra_args".
     - Logs stdout and stderr to a log file for debugging.
     - Uses a Singularity container for reproducibility and environment consistency.
     - Not included in rule all, as this is an optional evaluation step that may not be needed for all users.
    """
    input:
        assembly=get_assembly_for_type
    output:
        report=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast/report.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast",
        container=QC_CONTAINER,
        extra=lambda wildcards: config.get("assembly_eval", {}).get("extra_args", "")
    threads: config["threads"]
    log: "logs/assembly_eval_metaquast_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} quast.py \
            -o {params.outdir} \
            -t {threads} \
            {params.extra} \
            {input.assembly} > {log} 2>&1
        """


rule aggregate_best_mag_references:
    """
    Aggregates curated best MAG FASTAs into a single reference file for future MetaQUAST runs.
    Not included in rule all.
        - Takes a list of best MAG reference FASTA files from the config (best_mag_reference_inputs).
        - Concatenates them into a single FASTA file at the specified output location.
        - Logs any issues during concatenation to a log file for debugging.
        - If no best MAG reference inputs are configured, logs a message and exits with an error.
        - This rule allows users to easily create a custom reference set for evaluating their assemblies against high-quality MAGs.
        - The resulting reference file can be used in the metaquast_best_mag_reference_eval rule for more targeted assembly evaluation.
        - Users should ensure that the best MAG reference inputs are properly curated and relevant to their samples for meaningful evaluation results.
        - Not included in rule all since this is an optional step that depends on user-provided references.
        - The output reference file will be used in subsequent MetaQUAST evaluations if provided, but the workflow can still run without it if not needed.
        - This modular approach allows users to choose whether or not to include curated references in their assembly evaluation process based on their specific research questions and available data.
        - Overall, this rule provides flexibility for users to enhance their assembly evaluation with custom references while maintaining the ability to run the workflow without them if desired.
        - Users should check the logs for any issues during aggregation and ensure that the resulting reference file is correctly formatted and contains the expected sequences before using it in downstream analyses.
    """
    input:
        mags=get_best_mag_reference_inputs
    output:
        reference=f"{OUTPUT_DIR}/references/best_mags_reference.fasta"
    log: "logs/aggregate_best_mag_references.log"
    shell:
        """
        mkdir -p $(dirname {output.reference}) logs
        if [ -z "{input.mags}" ]; then
            echo "No best_mag_reference_inputs configured in config.yaml" > {log}
            exit 1
        fi
        cat {input.mags} > {output.reference} 2> {log}
        """


rule metaquast_best_mag_reference_eval:
    """
    Evaluates assemblies with MetaQUAST against curated best MAG references.
    Not included in rule all.
        - Uses the assembly file determined by get_assembly_for_type and the aggregated best MAG reference from the previous rule.
        - Outputs a report.tsv file in a separate output directory structure to distinguish it from the standard MetaQUAST evaluation.
        - Allows for extra arguments to be passed to MetaQUAST via the config file under "assembly_eval.extra_args".
        - Logs stdout and stderr to a log file for debugging.
        - Uses a Singularity container for reproducibility and environment consistency.
        - This rule provides an optional, more targeted assembly evaluation against high-quality MAG references, which can be particularly useful for assessing how well the assemblies capture known taxa of interest in the samples.
        - Users should ensure that the best MAG references are relevant to their samples for meaningful evaluation results.
        - The resulting report can be compared to the standard MetaQUAST evaluation to see how well the assemblies perform against both general and specific reference sets.
        - Not included in rule all since this is an optional step that depends on user-provided references and specific research questions.
        - Overall, this rule enhances the flexibility of the assembly evaluation process by allowing users to incorporate custom references while still maintaining the ability to run a standard evaluation without them if desired.
        - Users should check the logs for any issues during evaluation and ensure that the resulting report contains the expected metrics and comparisons based on the provided references.
        - This targeted evaluation can provide valuable insights into the quality of assemblies in relation to known high-quality MAGs, which can inform downstream analyses and interpretations of metagenomic data.
        - By comparing results from both evaluations, users can gain a more comprehensive understanding of their assemblies' performance across different reference sets and make informed decisions about which assemblies to use for further analyses based on their specific research goals and available data.
        - Ultimately, this rule allows users to tailor their assembly evaluation process to better suit their needs and enhance the interpretability of their results in the context of known microbial diversity represented by curated MAGs.
        - Users should consider including this step if they have access to relevant high-quality MAG references that can provide additional insights into their assemblies' performance beyond what is captured by a general reference set in a standard MetaQUAST evaluation.
    """
    input:
        assembly=get_assembly_for_type,
        reference=rules.aggregate_best_mag_references.output.reference
    output:
        report=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast_best_mags/report.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/metaquast_best_mags",
        container=QC_CONTAINER,
        extra=lambda wildcards: config.get("assembly_eval", {}).get("extra_args", "")
    threads: config["threads"]
    log: "logs/assembly_eval_metaquast_best_mags_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} metaquast.py \
            -r {input.reference} \
            -o {params.outdir} \
            -t {threads} \
            {params.extra} \
            {input.assembly} > {log} 2>&1
        """


rule map_short_reads_to_assembly:
    """
    Maps short reads to assemblies for mapping-based assembly performance stats.
        - Uses minimap2 for mapping and samtools for sorting and indexing the resulting BAM files.
        - Takes the assembly file determined by get_assembly_for_type and the short interleaved reads specified in the config.
        - Outputs a sorted BAM file and its index in a structured output directory.
        - Logs stdout and stderr to a log file for debugging.
        - Uses a Singularity container for reproducibility and environment consistency.
        - This mapping step allows users to assess how well their assemblies capture the original read data, which can provide insights into assembly quality and completeness beyond what is captured by reference-based evaluations like MetaQUAST.
        - Users can use the resulting BAM files for further analyses, such as coverage assessment, variant calling, or other downstream applications that require mapped reads.
        - Not included in rule all since this is an optional step that may not be needed for all users depending on their specific research questions and analysis goals.
        - Overall, this rule provides an additional layer of assembly evaluation based on read mapping, which can complement reference-based assessments and enhance the interpretability of assembly performance in the context of the original sequencing data.
        - Users should check the logs for any issues during mapping and ensure that the resulting BAM files are correctly formatted and contain the expected mappings before using them in downstream analyses.
        - This mapping-based evaluation can provide valuable insights into how well the assemblies represent the underlying read data, which can inform decisions about which assemblies to use for further analyses based on their performance in capturing the original sequencing information.
        - By combining results from both reference-based evaluations (like MetaQUAST) and mapping-based assessments, users can gain a more comprehensive understanding of their assemblies' quality and make informed decisions about their suitability for downstream applications based on their specific research goals and available data.
    """
    input:
        reads=lambda wildcards: config["input_reads"]["short_interleaved"].format(sample=wildcards.sample),
        assembly=get_assembly_for_type
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam.bai"
    params:
        outdir=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping",
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/assembly_eval_map_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        singularity exec {params.container} sh -c "set -euo pipefail; \
            export PATH=/opt/conda/envs/qc_env/bin:/opt/conda/envs/base_tools/bin:/opt/conda/bin:$PATH; \
            command -v minimap2 >/dev/null; \
            command -v samtools >/dev/null; \
            minimap2 -t {threads} -ax sr {input.assembly} {input.reads} | \
            samtools sort -@ {threads} -o {output.bam}; \
            samtools index -@ {threads} {output.bam}" > {log} 2>&1
        """


rule samtools_mapping_stats:
    """
    Produces samtools flagstat and idxstats from mapped BAMs.
        - Takes the sorted BAM file from the mapping step and generates flagstat and idxstats reports.
        - Outputs the reports in the same mapping directory structure for easy organization.
        - Logs stdout and stderr to a log file for debugging.
        - Uses the QC Singularity container for reproducibility and environment consistency.
        - This rule provides important mapping statistics that can be used to assess the quality of the read mapping to the assemblies, which can in turn provide insights into assembly quality and completeness.
        - Users can use the resulting flagstat and idxstats reports to evaluate the proportion of reads that mapped to the assemblies, the distribution of mappings across contigs, and other relevant metrics that can
    inform decisions about which assemblies to use for further analyses based on their performance in capturing the original read data.
        - Not included in rule all since this is an optional step that may not be needed for all users depending on their specific research questions and analysis goals.
        - Overall, this rule enhances the assembly evaluation process by providing detailed mapping statistics that complement reference-based assessments and help users understand how well their assemblies represent the underlying sequencing data.
        - Users should check the logs for any issues during the generation of mapping statistics and ensure that the resulting reports contain the expected metrics before using them for interpretation or downstream analyses.
    """
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/short_reads.sorted.bam.bai"
    output:
        flagstat=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/flagstat.txt",
        idxstats=f"{OUTPUT_DIR}/{{sample}}/assembly_eval/{{assembly_type}}/mapping/idxstats.txt"
    params:
        container=QC_CONTAINER
    threads: config["threads"]
    log: "logs/assembly_eval_samtools_{sample}_{assembly_type}.log"
    shell:
        """
        mkdir -p logs
        singularity exec {params.container} sh -c "set -euo pipefail; \
            export PATH=/opt/conda/envs/qc_env/bin:/opt/conda/envs/base_tools/bin:/opt/conda/bin:$PATH; \
            command -v samtools >/dev/null; \
            samtools flagstat -@ {threads} {input.bam} > {output.flagstat}; \
            samtools idxstats {input.bam} > {output.idxstats}" > {log} 2>&1
        """
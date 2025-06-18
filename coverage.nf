#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
Sequences in:${params.sequences}
Bedfile:${params.bedfile_exonwise}

"""

process coverage {
	publishDir "${params.output}", mode: 'copy', pattern: "*.counts.bed"
	input:
		tuple val (Sample), file(final_bam), file(final_bai)
		file (bedfile)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${final_bam} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${bedfile} -b ${Sample}.bed > ${Sample}.counts.bed
	"""
}

process coverview_run {
	executor="local"
	publishDir "${params.output}", mode: 'copy', pattern: "*.coverview_regions.csv"
	input:
		tuple val (Sample), file(final_bam), file(final_bai)
		file (bedfile)
	output:
		tuple val (Sample), file ("*.coverview_regions.csv")
	script:
	"""
	${params.coverview_path}/coverview -i ${final_bam} -b ${bedfile} -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}

workflow COVERAGE {
	take:
		final_bams_ch
		bedfile
	main:
	coverage(final_bams_ch, bedfile)
	coverview_run(final_bams_ch, bedfile)
}

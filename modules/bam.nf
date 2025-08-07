#!/usr/bin/nextflow

process TRIM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(read1), file(read2)
		path(ADAPTOR)
	output:
		tuple val(Sample), file("${Sample}_trim_R1.fastq"), file("${Sample}_trim_R2.fastq")
	script:
	"""
	fastp -i ${read1} -I ${read2} -o ${Sample}_trim_R1.fastq -O ${Sample}_trim_R2.fastq --adapter_fasta ${ADAPTOR} -w $task.cpus
	"""
}

process MAPBAM {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val(Sample), file(trim1), file(trim2)
		path (genome_dir)
		file (genome_fasta)
	output:
		tuple val(Sample), file ("${Sample}.bam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" \
	-M -t $task.cpus ${genome_dir}/${genome_fasta} ${trim1} ${trim2} | samtools sort -@ 8 -o ${Sample}.bam -
	"""
}

process MARK_DUPS {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(sortd_bam)
	output:
		tuple val(Sample), file("${Sample}_markdups.bam"), file("${Sample}_marked_dup_metrics.txt")
	script:
	"""
	picard MarkDuplicates \
		I=${sortd_bam} \
		O=${Sample}_markdups.bam \
		M=${Sample}_marked_dup_metrics.txt \
		TMP_DIR=${Sample}_tmp
	"""
}

process BQSR {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics)
		path (genome_dir)
		file (genome_fasta)
		path (SNPS)
		path (SNPS_index)
		path (INDELS)
		path (INDELS_index)
	output:
		tuple val(Sample), file("${Sample}_recal.table")
	script:
	"""
	gatk BaseRecalibrator \
		-I ${markdups_bam} \
		-R ${genome_dir}/${genome_fasta} \
		--known-sites ${SNPS} \
		--known-sites ${INDELS} \
		--bqsr-baq-gap-open-penalty 30.0 \
		-O ${Sample}_recal.table
	"""
}

process APPLY_BQSR {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_final.bam'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_final.bam.bai'
	tag "${Sample}"
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics), file(recal_table)
		path (genome_dir)
		path (genome_fasta)
	output:
		tuple val(Sample), file("${Sample}_final.bam"), file("${Sample}_final.bam.bai")
	script:
	"""
	gatk ApplyBQSR \
		-R ${genome_dir}/${genome_fasta} \
		-I ${markdups_bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${Sample}_final.bam

	mv ${Sample}_final.bai ${Sample}_final.bam.bai
	"""
}

process ALIGNMENT_METRICS {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_alignment_summary_metrics.txt'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
		path (genome_dir)
		path (genome_fasta)
	output:
		tuple val(Sample), file("${Sample}_alignment_summary_metrics.txt")
	script:
	"""
	picard CollectAlignmentSummaryMetrics \
		R=${genome_dir}/${genome_fasta} \
		I=${final_bam} \
		O=${Sample}_alignment_summary_metrics.txt
	"""
}

process INSERT_SIZE_METRICS {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.txt'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.pdf'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_insert_size_metrics.txt"), file("${Sample}_insert_size_metrics.pdf")
	script:
	"""
	picard CollectInsertSizeMetrics \
		I=${final_bam} \
		O=${Sample}_insert_size_metrics.txt \
		H=${Sample}_insert_size_metrics.pdf \
		HISTOGRAM_WIDTH=500 \
		TMP_DIR=${Sample}_tmp
	"""
}

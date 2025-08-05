#!/usr/bin/nextflow

process TRIM_RNA {
	label 'process_high'
	tag "${Sample}"
	input:
		tuple val(Sample), file(read1), file(read2)
		path (Adapt)
	output:
		tuple val(Sample), file("${Sample}_1P.fq.gz"), file("${Sample}_2P.fq.gz")
	script:
	"""
	trimmomatic PE \
	${read1} ${read2} \
	-threads $task.cpus -baseout ${Sample}.fq.gz ILLUMINACLIP:${Adapt}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	"""
}

process STAR {
	label 'process_high'
	tag "${Sample}"
	publishDir "${params.output}/${Sample}", mode : 'copy', pattern: '*.Aligned.sortedByCoord.out.bam'
	input:
		tuple val(Sample), file(ForwardRead), file(ReverseRead)
		file (STAR)
		file (STAR_GTF)
	output:
		tuple val(Sample), file("${Sample}.Aligned.sortedByCoord.out.bam")
	script:
	"""
	STAR --genomeDir ${STAR} --readFilesIn ${ForwardRead} ${ReverseRead} --runThreadN $task.cpus \
	--outFileNamePrefix ${Sample}. --sjdbGTFfile ${STAR_GTF} \
	--outSAMattrRGline ID:${Sample} 'SM:${Sample}' --twopassMode Basic \
	--chimOutType WithinBAM --chimSegmentMin 20 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 \
	--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
	"""
}

process INDEX {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}", mode : 'copy', pattern: '*.Aligned.sortedByCoord.out.bam.bai'
	input:
		tuple val (Sample), file(BamFile)
	output:
		tuple val (Sample), file("${Sample}.Aligned.sortedByCoord.out.bam.bai")
	script:
	"""
	samtools index ${BamFile} > ${Sample}.Aligned.sortedByCoord.out.bam.bai
	"""
}

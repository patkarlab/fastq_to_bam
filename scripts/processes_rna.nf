#!/usr/bin/nextflow

process Trim {
	label 'process_high'
	tag "${Sample}"
	container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
	input:
		tuple val(Sample), file(read1), file(read2)
		path (Adapt)
	output:
		tuple val(Sample), file("${Sample}_1P.fq.gz"), file("${Sample}_2P.fq.gz")
	script:
	"""
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz \
	-threads $task.cpus -baseout ${Sample}.fq.gz ILLUMINACLIP:${Adapt}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	"""
}

process AlignSTAR {
	label 'process_high'
	tag "${Sample}"
	publishDir "${params.output}/${Sample}", mode : 'copy', pattern: '*.Aligned.sortedByCoord.out.bam'
	container 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'
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

process Index {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}", mode : 'copy', pattern: '*.Aligned.sortedByCoord.out.bam.bai'
	container 'quay.io/biocontainers/samtools:1.15.1--h1170115_0'
	input:
		tuple val (Sample), file(BamFile)
	output:
		tuple val (Sample), file("${Sample}.Aligned.sortedByCoord.out.bam.bai")
	script:
	"""
	samtools index ${BamFile} > ${Sample}.Aligned.sortedByCoord.out.bam.bai
	"""
}

#!/usr/bin/nextflow
nextflow.enable.dsl=2

include { FASTQTOBAM } from './processes.nf'
include { COVERAGE } from './coverage.nf'

workflow FASTQ_TO_BAM{
	
	Channel
		.fromPath(params.input)
		.splitCsv(header: false)
		.map { row -> row[0] }
		.map { sample ->
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz")
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz")
			tuple(sample, r1, r2)
		}
		.set { samples_ch }

	FASTQTOBAM(samples_ch)
}

workflow FASTQ_TO_COV{
	
	Channel
		.fromPath(params.input)
		.splitCsv(header: false)
		.map { row -> row[0] }
		.map { sample ->
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz")
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz")
			tuple(sample, r1, r2)
		}
		.set { samples_ch }

	bedfile = file("${params.bedfile_exonwise}", checkIfExists: true)

	final_bams_ch = FASTQTOBAM(samples_ch)
	COVERAGE(final_bams_ch, bedfile)
}

workflow BAM_TO_COV {
	Channel
		.fromPath(params.input)
		.splitCsv(header: false)
		.map { row -> row[0] }
		.map { sample ->
			def bam = file("${params.sequences}/${sample}*.bam")
			def bai = file("${params.sequences}/${sample}*.bam.bai")
			tuple(sample, bam, bai)
		}
		.set { bams_ch }

	bedfile = file("${params.bedfile_exonwise}", checkIfExists: true)

	COVERAGE(bams_ch, bedfile)
}


workflow.onComplete {
		log.info ( workflow.success ? "\n\nDone! Output in the ${params.output} directory \n" : "Oops .. something went wrong" )
}

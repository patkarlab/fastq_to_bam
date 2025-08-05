#!/usr/bin/nextflow
nextflow.enable.dsl=2

// file paths
adapt = file("${params.adaptors}", checkIfExists: true)
star = file("${params.star_reference}", checkIfExists: true)
star_gtf = file("${params.star_gtf_path}", checkIfExists: true)

include { TRIM_RNA; STAR; INDEX; FEATURECOUNTS } from '../modules/bam_rna.nf'

workflow FASTQ_TO_BAM_RNA{	
	Channel
		.fromPath(params.input)
		.splitCsv(header: false)
		.map { sample ->
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { samples_ch }

	TRIM_RNA(samples_ch, adapt)
	STAR(TRIM_RNA.out, star, star_gtf)
	INDEX(STAR.out.star_bam)
	FEATURECOUNTS(STAR.out.star_bam, star_gtf)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the ${params.output} directory \n" : "Oops .. something went wrong" )
}
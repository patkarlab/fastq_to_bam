#!/usr/bin/nextflow
nextflow.enable.dsl=2

include { Trim } from './processes_rna.nf'
include { AlignSTAR } from './processes_rna.nf'
include { Index } from './processes_rna.nf'

workflow {
	
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

	adapt = file("${params.adaptors}", checkIfExists: true)
	star = file("${params.star_reference}")
	star_gtf = file("${params.star_gtf_path}")

	Trim(samples_ch, adapt)
	AlignSTAR(Trim.out, star, star_gtf)
	Index(AlignSTAR.out)
}

workflow.onComplete {
        log.info ( workflow.success ? "\n\nDone! Output in the ${params.output} directory \n" : "Oops .. something went wrong" )
}


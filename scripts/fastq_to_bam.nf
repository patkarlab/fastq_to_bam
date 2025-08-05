#!/usr/bin/nextflow
nextflow.enable.dsl=2

include { FASTQTOBAM } from './processes.nf'
include { ICHOR_CNA } from './ichorCNA.nf'
include { COVERAGE_WGS } from './coverage_lowpass_wgs.nf'

workflow {
	
	Channel
		.fromFilePairs("./sequences/*_S[0-9]*_R{1,2}_001.fastq.gz", flat: true)
		.set { samples_ch }

	final_bams_ch = FASTQTOBAM(samples_ch)
	ICHOR_CNA(final_bams_ch)
	COVERAGE_WGS(final_bams_ch)
}

workflow.onComplete {
        log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}


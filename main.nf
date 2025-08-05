#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_TO_BAM } from './workflows/fastq_bam'
include { FASTQ_TO_BAM_RNA } from './workflows/fastq_bam_rna'

//
// WORKFLOW: Run main nf-core/rnafusion analysis pipeline
//

workflow FASTQ_BAM {
	FASTQ_TO_BAM ()
}

workflow FASTQ_BAM_RNA {
	FASTQ_TO_BAM_RNA ()
}
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQ_TO_BAM } from './workflows/fastq_bam'
include { FASTQ_TO_BAM_RNA } from './workflows/fastq_bam_rna'

//
// WORKFLOW: Run main fastq to bam analysis pipeline
//

workflow FASTQ_BAM {
	FASTQ_TO_BAM ()
}

workflow FASTQ_BAM_RNA {
	FASTQ_TO_BAM_RNA ()
}
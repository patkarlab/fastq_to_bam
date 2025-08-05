#!/usr/bin/nextflow
nextflow.enable.dsl=2

// file paths
adapt = file("${params.adaptors}", checkIfExists: true )
seq_locatn = file("${params.sequences}", checkIfExists: true )
genome_loc = file("${params.genome}", checkIfExists: true)
genome_index = file("${params.genome_idx}", checkIfExists: true)
genome_dict = file("${params.genome_dict}", checkIfExists: true)
known_SNPs = file("${params.site1}", checkIfExists: true)
known_SNPs_index = file("${params.site1_idx}", checkIfExists: true)
known_INDELS = file("${params.site2}", checkIfExists: true)
known_INDELS_index = file("${params.site2_idx}", checkIfExists: true)

include { TRIM; MAPBAM; SORT; MARK_DUPS; BQSR; APPLY_BQSR; ALIGNMENT_METRICS; INSERT_SIZE_METRICS} from '../modules/bam.nf'
// include { COVERAGE } from './coverage.nf'

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
	TRIM(samples_ch, adapt)
	MAPBAM(TRIM.out, genome_loc)
	SORT(MAPBAM.out)
	MARK_DUPS(SORT.out)
	BQSR(MARK_DUPS.out, genome_loc, genome_index, genome_dict, known_SNPs, known_SNPs_index, known_INDELS, known_INDELS_index)
	APPLY_BQSR(MARK_DUPS.out.join(BQSR.out), genome_loc, genome_index, genome_dict)
	ALIGNMENT_METRICS(APPLY_BQSR.out, genome_loc, genome_index, genome_dict)
	INSERT_SIZE_METRICS(APPLY_BQSR.out)
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

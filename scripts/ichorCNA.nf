#!/usr/bin/nextflow

process RUN_ICHOR {
    tag "${Sample}"
    publishDir "${params.output}/${Sample}/" , mode: 'copy', pattern: '*_ichorCNA'

    input:
    tuple val(Sample), file(final_bam), file(final_bai)

    output:
    tuple val(Sample), path("${Sample}_ichorCNA")

    script:
    """
    ${params.run_ichorCNA} ${Sample} ${final_bam}

    """
}


workflow ICHOR_CNA {
    take:
        final_bams_ch
    main:
        RUN_ICHOR(final_bams_ch)
}

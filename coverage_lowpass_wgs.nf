#!/usr/bin/nextflow

process COVERAGE {
    tag "${Sample}"
    publishDir "${params.output}/${Sample}/" , mode: 'copy', pattern: '*txt'

    input:
    tuple val(Sample), file(final_bam), file(final_bai)

    output:
    tuple val(Sample), file("${Sample}_WGS_coverage_metrics.txt")

    script:
    """
    ${params.java_path}/java -jar ${params.picard_path} CollectWgsMetrics \
    I=${final_bam} \
    R=${params.genome} \
    O=${Sample}_WGS_coverage_metrics.txt
    """
}


workflow COVERAGE_WGS {
    take:
        final_bams_ch
    main:
        COVERAGE(final_bams_ch)
}

#!/bin/bash

#nextflow -c /home/diagnostics/pipelines/fastq_to_bam_paired_snakemake/nextflow.config run fastq_to_bam.nf --output Final_Output -resume -bg -with-dag

nextflow -c /home/diagnostics/pipelines/fastq_to_bam_paired_snakemake/nextflow.config run fastq_to_bam.nf --genome /home/reference_genomes/hg38_broad/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --site1 /home/reference_genomes/dbSNPGATK_hg38/Homo_sapiens_assembly38.dbsnp138.vcf --site2 /home/reference_genomes/dbSNPGATK_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --run_ichorCNA /home/diagnostics/pipelines/fastq_to_bam_paired_snakemake/scripts/run_ichorCNA_v4.sh --output Final_Output -resume -bg

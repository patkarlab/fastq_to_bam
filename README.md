# fastq_to_bam

This is a nextflow pipeline to convert paired-end fastqs to bams, following GATK best practices. It is based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake

For running this pipeline, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- sequences = path to the folder containing fastq files

- genome = path to the genomic fasta file
    - eg. `/home/reference_genomes/hg19_broad/hg19_all.fasta`
    - This folder should also contain .fai and .dict files for the genome 
		- eg. `/home/reference_genomes/hg19_broad/hg19_all.fasta.fai`
        - eg. `/home/reference_genomes/hg19_broad/hg19_all.fasta.dict`

- site1 = path to the vcf file containing known polymorphisms. (dbsnp_138)
    - eg. `/home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.vcf`
    - This folder should also contain .idx file for the vcf.
        - eg. `/home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.vcf.idx`

- site2 = path to the vcf file containing known indels. (Mills_and_1000G_gold_standard.indels)
    - eg. `/home/reference_genomes/dbSNPGATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf`
    - This folder should also contain .idx file for the vcf.
    	- eg. `/home/reference_genomes/dbSNPGATK Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx`

- adaptors = path to the Fasta file of adapter sequences for trimming 
    - eg. `./scripts/TruSeq2-PE.fa`

- star_reference = path to the STAR reference folder
     - eg. `home/diagnostics/pipelines/nf-core/rnafusion/references/star`

- star_gtf_path = path to the GTF file for STAR
     - eg. `/home/diagnostics/pipelines/nf-core/rnafusion/references/ensembl/Homo_sapiens.GRCh38.102.gtf`

## Usage:

### For DNA:

```
nextflow run main.nf -entry FASTQ_BAM -profile docker -resume -bg
```

### For RNA:

```
nextflow run main.nf -entry FASTQ_BAM_RNA -profile docker -resume -bg
```

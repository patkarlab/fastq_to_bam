# fastq_to_bam

This is a nextflow pipeline to convert paired-end fastqs to bams, following GATK best practices. It is based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake

For running this pipeline, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- genome = Genomic fasta file.
- site1 = known_polymorphic_sites 1 (dbsnp_138)
- site2 = known_polymorphic_sites 2 (Mills_and_1000G_gold_standard.indels)
- adaptors = Fasta file of adapter sequences for trimming 
- java_path = directory containing the java executable
- samtools = samtools executable path
- picard_path = path to the picard.jar file
- gatk = gatk executable path
- fastp = fastp executable path
- star_reference = STAR reference folder
- star_gtf_path = GTF file for STAR
- bwa - bwa executable path


## Usage:

### For DNA:

```
get_bam -i <fastq_location> -s <samples_file>

    -i <fastq_location>        Complete path to the directory containing FASTQ files
    -s <samples_file>          List of the sample names
    -o <output>                Output folder. If not specified, output will be written to the BAM folder at the fastq location

Options:
    -h, --help                 Display this help message and exit
```

### For RNA:

```
get_bam_rna -i <fastq_location> -s <samples_file>

   -i <fastq_location>        Complete path to the directory containing FASTQ files
   -s <samples_file>          List of the sample names
   -o <output>                Output folder. If not specified, output will be written to the BAM folder at the fastq location

Options:
   -h, --help                 Display this help message and exit
```

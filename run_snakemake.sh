#!/usr/bin/bash


snakemake -p -s fastq_to_bam_paired.snakefile --until get_insert_size_metrics --cores 8 --dag dot | dot -Tsvg > dag.svg

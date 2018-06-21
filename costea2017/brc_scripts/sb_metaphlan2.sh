#!/bin/bash
#SBATCH -o sb_metaphlan2_%j.out
#SBATCH -c 8

accession=$1
reads_path=/home/mrmclare/data/costea2017/reads
out_path=/home/mrmclare/data/costea2017/metaphlan2

/home/mrmclare/applications/biobakery-metaphlan2-5175d8783801/metaphlan2.py \
    --bowtie2out $out_path/${accession}_bowtie2.bz2 \
    --nproc 10 --input_type fastq \
    $reads_path/${accession}_1.fastq.gz,$reads_path/${accession}_2.fastq.gz \
    $out_path/${accession}_profiled_metagenome.tsv

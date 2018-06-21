#!/bin/bash
#SBATCH -o sb_metaphlan2_%j.out
#SBATCH -c 4

accession=$1
reads_path=/home/mrmclare/data/costea2017/reads
out_path=/home/mrmclare/data/costea2017/metaphlan2

/home/mrmclare/applications/biobakery-metaphlan2-5175d8783801/metaphlan2.py \
    -t clade_profiles \
    --nproc 4 --input_type bowtie2out \
    $out_path/${accession}_bowtie2.bz2 \
    $out_path/${accession}_clade_profiles.txt

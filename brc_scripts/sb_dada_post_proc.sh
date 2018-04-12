#!/bin/bash
#SBATCH -o sb_dada_post_proc_%j.out
#SBATCH -c 10
Rscript ~/metagenomics_calibration/mbqc/dada2_pipeline/brc_scripts/dada_post_proc.R

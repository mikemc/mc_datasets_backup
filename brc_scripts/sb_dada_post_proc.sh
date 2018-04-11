#!/bin/bash
#SBATCH -o sb_dada_post_proc_%j.out
#SBATCH -c 10
Rscript ~/metagenomics_calibration/mbqc/brc_scripts/dada_post_proc.R

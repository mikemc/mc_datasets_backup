#!/bin/bash
#SBATCH -o sb_assign_taxonomy_%j.out
#SBATCH -c 10
Rscript ~/metagenomics_calibration/mbqc/dada2_pipeline/brc_scripts/assign_taxonomy.R

#!/bin/bash
#SBATCH -o sb_dada_pipeline_%j.out
#SBATCH -c 10
Rscript ~/metagenomics_calibration/mbqc/brc_scripts/dada_pipeline.R $1

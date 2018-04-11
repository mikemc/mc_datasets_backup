#!/bin/bash
#SBATCH -o sb_assign_species_%j.out
#SBATCH -c 10
Rscript ~/metagenomics_calibration/mbqc/brc_scripts/assign_species.R

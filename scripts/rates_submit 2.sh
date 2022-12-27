#!/bin/bash
#SBATCH -c 4
#SBATCH --time 10:00:00
module load R/4.2.0-foss-2021b
Rscript extract_validation_sims_rates_large_indiv.Rscript $1 

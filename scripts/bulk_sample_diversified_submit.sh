#!/bin/sh
#SBATCH -c 6
#SBATCH --time 10:00:00
module load R/4.2.0-foss-2021b
Rscript bulk_sample_simulated_tumor_large_indiv_diversified.Rscript $1 

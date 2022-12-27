#!/bin/bash
#SBATCH --time 10:00:00
module load R/4.2.0-foss-2021b
Rscript get_tree_stats_indiv.Rscript $1 

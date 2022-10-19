#!/bin/bash
#SBATCH --time 10:00:00
module load R/4.2.0-foss-2021b
Rscript get_edge_dist_vs_divisions_indiv.Rscript  $1 

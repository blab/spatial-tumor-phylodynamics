#!/bin/bash
#SBATCH --time 10:00:00
module load R/4.2.0-foss-2021b
echo $1
echo $2
Rscript set_up_simulated_tumors_xmls_large_indiv_random.Rscript $1 $2

#!/bin/bash
#SBATCH -c 4
python3 -u run_eden_simulations_indiv_lowmem_small.py $1

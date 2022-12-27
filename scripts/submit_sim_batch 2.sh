#!/bin/sh
#SBATCH -c 10
echo "running"
python3 $1

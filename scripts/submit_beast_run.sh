#!/bin/bash
#SBATCH -c 2
module load Java/15.0.1
java -jar Nab7.jar $1

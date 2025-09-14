#!/bin/bash

#SBATCH --job-name=symmetric
#SBATCH --time=10:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=logs/symmetric.error
#SBATCH --output=logs/symmetric.log
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript test_symmetric.R
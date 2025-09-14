#!/bin/bash

#SBATCH --job-name=platinum_bootstrap
#SBATCH --time=50:00:00       			
#SBATCH --ntasks=3
#SBATCH --mem=3G       				
#SBATCH --error=logs/platinum_bootstrap.error
#SBATCH --output=logs/platinum_bootstrap.log
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript platinum_bootstrap.R
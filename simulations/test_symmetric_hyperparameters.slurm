#!/bin/bash

#SBATCH --job-name=symmetric_hyperparameters
#SBATCH --time=40:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=logs/symmetric_hyperparameters.error
#SBATCH --output=logs/symmetric_hyperparameters.log
#SBATCH -p shared

module load R/4.2.2-fasrc01
Rscript test_symmetric_hyperparameters.R
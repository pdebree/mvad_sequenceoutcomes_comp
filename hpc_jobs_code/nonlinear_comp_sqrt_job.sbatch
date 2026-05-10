#!/bin/bash
#t
#SBATCH --job-name=sqnonlinear_comp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=8GB
#SBATCH --account=torch_pr_774_general
#SBATCH --output=nonlinear_comp_sqrt_%A_%a.out
#SBATCH --error=nonlinear_comp_sqrt_%A_%a.err

module purge
module load r/4.5.1

R CMD BATCH --no-save --no-restore mvad_nonlinear_comp_sqrt.R

#!/bin/bash
#t
#SBATCH --job-name=nonlinear_comp_cv
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=120:00:00
#SBATCH --mem=128GB
#SBATCH --account=torch_pr_774_general
#SBATCH --output=nonlinear_comp_%A_%a.out
#SBATCH --error=nonlinear_comp_%A_%a.err

module purge
module load r/4.5.1

R CMD BATCH --no-save --no-restore mvad_nonlinear_comp_cv.R

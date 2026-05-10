#!/bin/bash
#t
#SBATCH --job-name=linear_comp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=5GB
#SBATCH --account=torch_pr_774_general
#SBATCH --output=linear_comp_%A_%a.out
#SBATCH --error=linear_comp_%A_%a.err

module purge
module load r/4.5.1

R CMD BATCH --no-save --no-restore mvad_linear_comp.R

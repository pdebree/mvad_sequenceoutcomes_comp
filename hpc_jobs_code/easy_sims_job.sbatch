#!/bin/bash
#t
#SBATCH --job-name=easy_sims
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=120:00:00
#SBATCH --mem=128GB
#SBATCH --account=torch_pr_774_general
#SBATCH --output=easy_sims_%A_%a.out
#SBATCH --error=easy_sims_%A_%a.err

module purge
module load r/4.5.1

R CMD BATCH --no-save --no-restore easy_sim_comp.R

#!/bin/bash
#t
#SBATCH --job-name=dar_med_sims
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=120:00:00
#SBATCH --mem=128GB
#SBATCH --account=torch_pr_774_general
#SBATCH --output=dar_med_sims_%A_%a.out
#SBATCH --error=dar_med_sims_%A_%a.err

module purge
module load r/4.5.1

R CMD BATCH --no-save --no-restore dar_med_sim_comp.R

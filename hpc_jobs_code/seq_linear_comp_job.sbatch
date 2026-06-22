#!/bin/bash
#t
#SBATCH --job-name=seq_linear
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=64GB
#SBATCH --account=torch_pr_774_general
#SBATCH --output=linear_seq_%A_%a.out
#SBATCH --error=linear_seq_%A_%a.err

module purge
module load r/4.5.1

R CMD BATCH --no-save --no-restore mvad_linear_seqmets.R

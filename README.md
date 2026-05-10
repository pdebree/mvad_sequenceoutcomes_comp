# Sequence Outcomes Files

Files: 
- seqout_utils.R - Contains functions used in Linear, Non-Linear and Simulated Data Competitions. 
- comp_plots.R - Code for making plots based on Linear and Non-Linear Competitions
- mvad_demos_only.R - Code to evaluate the predictive power of demographics only
- mvad_linear_comp.R - Contains full linear competition (including sequence metrics work) 
- mvad_linear_comp_cv.R - Contains full linear competition (not including sequence metrics) to be run on HPC with 20 cvs
- mvad_nonlinear_comp.R - Contains full nonlinear competition (including sequence metrics work) 
- mvad_nonlinear_comp_cv.R - Contains full nonlinear competition (not including sequence metrics) to be run on HPC with 20 cvs
- mvad_varex_comparisons.R - Code for looking at decompositions without prediction (for comparison against choice from competition)
- easy_sim_comp.R - Full Competition run on "easy" simulated data 
- med_sim_comp.R - Full Competition run on "medium" simulated data 
- dar_med_sim_comp.R - Full Competition run on "DAR medium" simulated data 
- hard_sim_comp.R - Full Competition run on "hard" simulated data 
- very_hard_sim_comp.R - Full Competition run on "very hard" simulated data 
- sims_output_checker.qmd - Code for checking the outputs of simulations and making additional plots. 


Directories:
- plots - Directory containing plots created for analysis
- hpc_jobs_code - sbatch scripts for HPC runs 
- cv_outputs - contains files for cross validation of linear/non-linear competition outputs
- sims_outputs/conv_sims - contains files for sim study outputs




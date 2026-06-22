# Code to re-write simulations plots 

# Pippi de Bree
# gcd2056@nyu.edu
# 2026-06-22


library(tidyverse)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)


# Read in MSES
easy_mses <- readRDS("sims_output/sims_final_run/easy_seqs_full_mses.rds")
med_mses <- readRDS("sims_output/sims_final_run/med_seqs_full_mses.rds")
dar_med_mses <- readRDS("sims_output/sims_final_run/dar_med_seqs_full_mses.rds")
hard_mses <- readRDS("sims_output/sims_final_run/hard_seqs_full_mses.rds")
very_hard_mses <- readRDS("sims_output/sims_final_run/very_hard_seqs_full_mses.rds")


method_names <- c("CFDA", "Counts", "OM-Trate (Hard)","OM-SLOG (Hard)",
  "LCS (Hard)", "OM-Trate (Soft)","OM-SLOG (Soft)", "LCS (Soft)")


make_mses_plot <- function(mse.sims, sim_difficulty) {


  plot_data <- map_df(names(mse.sims), function(name) {
    matrix_data <- mse.sims[[name]]
    row_means <- sqrt(rowMeans(matrix_data))

    tibble(
      Method = factor(method_names, levels = method_names),
      Mean_MSE = row_means,
      Simulation_Type = name
    )
  })



  p <- ggplot(plot_data, aes(x = Method, y = Mean_MSE, color = Simulation_Type, group = Simulation_Type)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_minimal(base_size = 14) + 
    labs(
      title = paste0("Mean RMSE per Method - ", sim_difficulty, " Difficulty"),
      x = "Method",
      y = "Average RMSE (of Best Number of Components)",
      color = "Simulation Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


  return(p)
}

easy_p <- make_mses_plot(easy_mses, "Easy")
pdf("plots/easy_sim_mse_plot.pdf" ,width=9,height=7)
easy_p
dev.off()


med_p <- make_mses_plot(med_mses, "Medium")
pdf("plots/med_sim_mse_plot.pdf" ,width=9,height=7)
med_p
dev.off()

dar_med_p <- make_mses_plot(dar_med_mses, "Medium (DAR)")
pdf("plots/dar_med_sim_mse_plot.pdf" ,width=9,height=7)
dar_med_p 
dev.off()

hard_p <- make_mses_plot(hard_mses, "Hard")
pdf("plots/hard_sim_mse_plot.pdf" ,width=9,height=7)
hard_p
dev.off()


very_hard_p <- make_mses_plot(very_hard_mses, "Hard")
pdf("plots/very_hard_sim_mse_plot.pdf" ,width=9,height=7)
very_hard_p
dev.off()


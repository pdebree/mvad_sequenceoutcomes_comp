# Code to look at the difference in performance with different number of priors
# For soft clustering, try with 50, 100 and 200 priors

# Pippi de Bree
# gcd2056@nyu.edu
# 2026-06-18


library(TraMineR)
library(TraMineRextras)
library(cluster)
library(tidyverse)
library(fpc)
library(NbClust)
library(cfda)
library(ranger)
library(tuneRanger)
library(gt)
library(dplyr)
library(factoextra)
library(vivid)
source("seqout_utils.R")


nonlinear_comp <- readRDS("cv_outputs/CV_NonLinearComp.rds")
nonlinear_comp_sqrt <- readRDS("cv_outputs/sensitivity_analysis_runs/CV_NonLinearComp_sqrt.rds")


methods <- c("om_trate_hard", "om_trate_soft", "om_slog_hard", 
             "om_slog_soft", "lcs_hard", "lcs_soft", "windows", "harm")

method_labels <- c(
  "om_trate_hard" = "OM T-Rate (Hard)", 
  "om_trate_soft" = "OM T-Rate (Soft)", 
  "om_slog_hard"  = "OM INDELSLOG (Hard)", 
  "om_slog_soft"  = "OM INDELSLOG (Soft)",
  "lcs_hard"      = "LCS (Hard)", 
  "lcs_soft"      = "LCS (Soft)", 
  "windows"       = "Counts", 
  "harm"          = "CFDA"
)


get_nonlinear_metrics <- function(data_list, method_name) {
  
  # Extract the list of 20 CV 3d arrays (fold, comps, mtry)
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)
  
  # within every cross validation - find mse average by mtry 
  avg_cv_mse <- apply(combined_mse, c(2, 3, 4), mean, na.rm = TRUE)

  # essentially the location maps to the mtry values
  best_mtry_cv <- apply(avg_cv_mse, c(1, 3), safe_which_min)
  modal_mtrys <- apply(best_mtry_cv, 1, safe_mode)

  # the actual min corresponds to this so we will get mses that correspond to the best
  best_avg_mse_comps_cv <- apply(avg_cv_mse, c(1, 3), safe_min)
  # now we can average these mses (picked by mtry) to find mse by comp
  rmse_comps <- sqrt(apply(best_avg_mse_comps_cv, 1, mean, na.rm=TRUE))

  n_converged_all <- 100 - colSums(apply(combined_mse, c(2,3), is.na))
  n_converged <- n_converged_all[modal_mtrys]

  return(list(rmses = rmse_comps, modal_mtry = modal_mtrys, mtry_dis = best_mtry_cv, method_name = method_name, n_converged=n_converged))

}




# nonlinear with up n-1 mtry
nonlin_maxmin <- map_dfr(names(method_labels), function(m) {
  res <- get_nonlinear_metrics(nonlinear_comp, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmses,
    method = method_labels[[m]], 
    modal_mtry = res$modal_mtry,
    mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis))), 
    n_converged = res$n_converged
  )
})
nonlin_maxmin$mtry_max <- "maxmin"



# nonlinear with up n-1 mtry
nonlin_sqrt <- map_dfr(names(method_labels), function(m) {
  res <- get_nonlinear_metrics(nonlinear_comp_sqrt, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmses,
    method = method_labels[[m]], 
    modal_mtry = res$modal_mtry,
    mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis))), 
    n_converged = res$n_converged
  )
})
nonlin_sqrt$mtry_max <- "sqrt"



nonlinear_mtry_eval <- rbind(nonlin_maxmin, nonlin_sqrt)



ggplot(nonlinear_mtry_eval, aes(x = comps, y = rmse, color = factor(mtry_max), group = factor(mtry_max))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  scale_x_continuous(breaks = seq(2, 25, by = 2)) +
  scale_color_discrete(drop = FALSE) + 
  facet_wrap(~ method) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.background = element_rect(fill = "white", color = "grey80"),
    strip.text = element_text(face = "bold", size = 11)
  ) + 
  labs(
    title = "Random Forest Performance by Max Mtry (Non Linear)",
    x = "Number of Components / Clusters",
    y = "RMSE (Squared Average MSE Across CVs)",
    color = "Max Mtry Type"
  )



ggplot(nonlinear_mtry_eval, aes(x = comps, y = modal_mtry, color = factor(mtry_max), group = factor(mtry_max))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, alpha=0.5) + 
  scale_x_continuous(breaks = seq(2, 25, by = 2)) +
  scale_color_discrete(drop = FALSE) +  # Forces 1.25 onto the legend scale if it exists
  facet_wrap(~ method) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.background = element_rect(fill = "white", color = "grey80"),
    strip.text = element_text(face = "bold", size = 11)
  ) + 
  labs(
    title = "Random Forest Performance by Max Mtry (Non Linear)",
    x = "Number of Components / Clusters",
    y = "Modal Mtry",
    color = "Max Mtry Type"
  )

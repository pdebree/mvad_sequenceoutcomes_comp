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


linear_comp_100 <- readRDS("cv_outputs/CV_LinearComp.rds")
nonlinear_comp_100 <- readRDS("cv_outputs/CV_NonLinearComp.rds")

linear_comp_50 <- readRDS("cv_outputs/sensitivity_analysis_runs/CV_LinearComp_50prior.rds")
nonlinear_comp_50 <- readRDS("cv_outputs/sensitivity_analysis_runs/CV_NonLinearComp_prior50.rds")

linear_comp_200 <- readRDS("cv_outputs/sensitivity_analysis_runs/CV_LinearComp_200prior.rds")
nonlinear_comp_200 <- readRDS("cv_outputs/sensitivity_analysis_runs/CV_NonLinearComp_prior200.rds")

# MSE Plot
get_method_metrics <- function(data_list, method_name, method_labs) {

  # pull out cross validated mses for given method
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)

  # find means across the folds (of all cross )
  mean_mse <- apply(combined_mse, 2, mean, na.rm = TRUE)


  df <- data.frame(
      Method = method_name, 
      Index = 1:length(mean_mse),
      RMSE   = sqrt(mean_mse), 
      stringsAsFactors = FALSE, 
      n_converged = 100 - colSums(apply(combined_mse, 2, is.na))
    )


  df <- df |> mutate(Method = recode(Method, !!!method_labs))

  return(df)
}


soft_methods <- c("om_trate_soft", "om_slog_soft", "lcs_soft")

soft_method_labels <- c(
  "om_trate_soft" = "OM T-Rate (Soft)", 
  "om_slog_soft"  = "OM INDELSLOG (Soft)",
  "lcs_soft"      = "LCS (Soft)"
)


# Nonlinear from regular run (1.5)
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





# linear with 100 priors
soft_lin_100 <- do.call(rbind, lapply(soft_methods, function(m) get_method_metrics(linear_comp_100, m, soft_method_labels)))
soft_lin_100$n_priors <- 100


# nonlinear with 100 priors
soft_nonlin_100 <- map_dfr(names(soft_method_labels), function(m) {
  res <- get_nonlinear_metrics(nonlinear_comp_100, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmses,
    method = soft_method_labels[[m]], 
    modal_mtry = res$modal_mtry,
    mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis))), 
    n_converged = res$n_converged
  )
})
soft_nonlin_100$n_priors <- 100



# linear with 50 priors
soft_lin_50 <- do.call(rbind, lapply(soft_methods, function(m) get_method_metrics(linear_comp_50, m, soft_method_labels)))
soft_lin_50$n_priors <- 50


# nonlinear with 1.25
soft_nonlin_50 <- map_dfr(names(soft_method_labels), function(m) {
  res <- get_nonlinear_metrics(nonlinear_comp_50, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmses,
    method = soft_method_labels[[m]], 
    modal_mtry = res$modal_mtry,
    mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis))), 
    n_converged = res$n_converged
  )
})
soft_nonlin_50$n_priors <- 50




# linear with 200
soft_lin_200 <- do.call(rbind, lapply(soft_methods, function(m) get_method_metrics(linear_comp_200, m, soft_method_labels)))
soft_lin_200$n_priors <- 200


# nonlinear with 200
soft_nonlin_200 <- map_dfr(names(soft_method_labels), function(m) {
  res <- get_nonlinear_metrics(nonlinear_comp_200, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmses,
    method = soft_method_labels[[m]], 
    modal_mtry = res$modal_mtry,
    mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis))), 
    n_converged = res$n_converged
  )
})
soft_nonlin_200$n_priors <- 200



# combine linear and nonlinear 
linear_npriors <- rbind(soft_lin_100, soft_lin_50, soft_lin_200)
nonlinear_npriors <- rbind(soft_nonlin_100, soft_nonlin_50, soft_nonlin_200)


pdf("plots/LinearPriorsEval.pdf",width=8,height=6)
ggplot(linear_npriors, aes(x = Index, y = RMSE, color = factor(n_priors), group = n_priors)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  scale_x_continuous(breaks = seq(2, 25, by = 2)) +
  facet_wrap(~ Method) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.background = element_rect(fill = "white", color = "grey80"),
    strip.text = element_text(face = "bold", size = 11)
  ) + 
  labs(
    title = "Random Forest Performance by Number of Priors (Linear)",
    x = "Number of Components / Clusters",
    y = "RMSE (Squared Average MSE Across CVs)",
    color = "Number of Priors"
  )
dev.off()

pdf("plots/NonLinearPriorsEval.pdf",width=8,height=6)
ggplot(nonlinear_npriors, aes(x = comps, y = rmse, color = factor(n_priors), group = factor(n_priors))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
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
    title = "Random Forest Performance by Number of Priors (Non Linear)",
    x = "Number of Components / Clusters",
    y = "RMSE (Squared Average MSE Across CVs)",
    color = "Number of Priors"
  )
dev.off()
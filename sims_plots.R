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


easy_results <- readRDS("cv_outputs/soft_6_sims/easy_seqs_results_10_175.rds")
med_results <- readRDS("cv_outputs/soft_6_sims/dar_med_seqs_results_10_175.rds")
hard_results <- readRDS("cv_outputs/soft_6_sims/hard_seqs_results_10_125.rds")


method_names <- c("CFDA", "Counts", "OM-Trate (Hard)","OM-SLOG (Hard)",
  "LCS (Hard)", "OM-Trate (Soft)","OM-SLOG (Soft)", "LCS (Soft)")


mses_plot_data <- function(mse.sims, diff_type) {
  plot_data <- map_df(names(mse.sims), function(name) {
    matrix_data <- mse.sims[[name]]
    row_means <- rowMeans(matrix_data, na.rm=TRUE)

    tibble(
      Method = factor(method_names, levels = method_names),
      Mean_MSE = row_means,
      Simulation_Type = name
    )
  })

  simulation_types <- list(
    "indep"="Indep.", 
    "semi"="Semi-I", 
     "semi_concord"="Semi-II",
      "concord"="Concordant")

  
  method_types <- list("CFDA" = "CFDA", 
  "Counts"="Counts-Based", 
  "OM-SLOG (Hard)" = "Hard Clust.",  
  "OM-SLOG (Soft)" = "Soft Clust.")
  
  
  
plot_data <- plot_data |> mutate(Simulation_Type = recode(Simulation_Type, !!!simulation_types)) |> 
  filter(Simulation_Type %in% simulation_types, Method %in% names(method_types)) |> mutate(Method = recode(Method, !!!method_types))
plot_data$Difficulty <- diff_type
  
return(plot_data)

}


easy_plot_data <- mses_plot_data(easy_results$mse, "Easy")
med_plot_data <- mses_plot_data(med_results$mse, "Medium")
hard_plot_data <- mses_plot_data(hard_results$mse, "Hard")

sims_plot_data <- rbind(easy_plot_data, med_plot_data, hard_plot_data)
sims_plot_data$Difficulty <- factor(sims_plot_data$Difficulty, levels = c("Easy", "Medium", "Hard"))
sims_plot_data$Simulation_Type <- factor(sims_plot_data$Simulation_Type, levels = c("Indep.","Semi-I", "Semi-II","Concordant"))

pdf("plots/SimulationResults.pdf" ,width=8,height=6)
ggplot(sims_plot_data, aes(x = Method, y = Mean_MSE, color = Simulation_Type, group = Simulation_Type)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_minimal(base_size = 14) + 
    labs(
      title = "RMSE by Feature Extraction Method - Simulated Data",
      subtitle = "DGP: Cluster Detection Level",
      x = "Method",
      y = "RMSE",
      color = "DGP: Link"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.subtitle = element_text(size = 12) 
  ) + facet_wrap(~ Difficulty)
dev.off()



# Code to check the convergences and number of condition numbers met 
hard_sims_train_conv <- hard_results$train_conv_tracker
med_sims_train_conv <- med_results$train_conv_tracker
easy_sims_train_conv <- easy_results$train_conv_tracker

hard_sims_test_conv <- hard_results$test_conv_tracker
med_sims_test_conv <- med_results$test_conv_tracker
easy_sims_test_conv <- easy_results$test_conv_tracker


non_convergence_train_summary <- function(train_conv_tracker) {  
  map_df(names(train_conv_tracker), function(sim_difficulty) {
  
  sets_list <- train_conv_tracker[[sim_difficulty]]

  map_df(seq_along(sets_list), function(set_idx) {
    methods_list <- sets_list[[set_idx]]
    
    map_df(names(methods_list), function(method_name) {

      if (method_name != "om_trate") {

      conv_matrix <- methods_list[[method_name]]

      conv_counts <- colSums(conv_matrix[,2:4,1])
      safe_condition_no <- colSums(conv_matrix[,2:4,2] < 100)
      mean_condition_no <- colMeans(conv_matrix[,2:4,2])
        
      tibble(
        Sim_Type = sim_difficulty,
        Method = method_name,
        Clusters = 2:4,
        Convergences = conv_counts,
        safe_condition_no = safe_condition_no, 
        mean_condition_no = mean_condition_no
      )
      }
      })
  })
}) %>%

  group_by(Sim_Type, Method, Clusters) %>%
  summarise(
    Total_Convergences = sum(Convergences), 
    Mean_Convergences = mean(Convergences),
    Total_Safe_Condition_No = sum(safe_condition_no), 
    Mean_Safe_Condition_No = mean(mean_condition_no),
  
  .groups = "drop") %>%

  mutate(
    Method = case_when(
      Method == "om_trate" ~ "OM-Trate (Soft)",
      Method == "om_slog"  ~ "OM-SLOG (Soft)",
      Method == "lcs"      ~ "LCS (Soft)",
      TRUE ~ Method
    )
  )  |> filter(Method == "OM-SLOG (Soft)", Clusters > 1, Sim_Type %in% 
    c("semi", "indep", "concord", "semi_concord"))

}



non_convergence_test_summary <- function(test_conv_tracker) {
  
  map_df(names(test_conv_tracker), function(sim_difficulty) {
  
  sets_list <- test_conv_tracker[[sim_difficulty]]
  
  map_df(seq_along(sets_list), function(set_idx) {
    methods_list <- sets_list[[set_idx]]
    
    map_df(names(methods_list), function(method_name) {      

      tibble(
        Sim_Type = sim_difficulty,
        Method = method_name,
        Train_Result_Available = methods_list[[method_name]][1] == 1, 
        Test_Converged = methods_list[[method_name]][2] == 1,
        Test_Condition_No = methods_list[[method_name]][3],
        Test_Safe_Condition_No = methods_list[[method_name]][3] < 100
      )
    })
  })
}) %>%
  group_by(Sim_Type, Method) %>%
  summarise(valid_train_res = sum(Train_Result_Available), 
    test_converged = sum(Test_Converged),
    test_safe_condition_no = sum(Test_Safe_Condition_No),
  .groups = "drop") %>%
  mutate(
    Method = case_when(
      Method == "om_trate" ~ "OM-Trate (Soft)",
      Method == "om_slog"  ~ "OM-SLOG (Soft)",
      Method == "lcs"      ~ "LCS (Soft)",
      TRUE ~ Method
    )
  ) |> filter(Method == "OM-SLOG (Soft)", Sim_Type %in% c("semi", "indep", "concord", "semi_concord"))
}


easy_train_non_conv_sum <- non_convergence_train_summary(easy_sims_train_conv)
easy_train_non_conv_sum$diff <- "Easy"
med_train_non_conv_sum <- non_convergence_train_summary(med_sims_train_conv)
med_train_non_conv_sum$diff <- "Medium"
hard_train_non_conv_sum <- non_convergence_train_summary(hard_sims_train_conv)
hard_train_non_conv_sum$diff <- "Hard"

train_conv_sum <- rbind(easy_train_non_conv_sum, med_train_non_conv_sum, hard_train_non_conv_sum)
print(train_conv_sum, n=Inf)

# Testing Conv Summary 
easy_test_non_conv_sum <- non_convergence_test_summary(easy_sims_test_conv)
easy_test_non_conv_sum$diff <- "Easy"
med_test_non_conv_sum <- non_convergence_test_summary(med_sims_test_conv)
med_test_non_conv_sum$diff <- "Medium"
hard_test_non_conv_sum <- non_convergence_test_summary(hard_sims_test_conv)
hard_test_non_conv_sum$diff <- "Hard"

test_conv_sum <- rbind(easy_test_non_conv_sum, med_test_non_conv_sum, hard_test_non_conv_sum)
print(test_conv_sum, n=Inf)

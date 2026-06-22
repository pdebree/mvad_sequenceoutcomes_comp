# Code for MVAD Sequence Work - Simulated Data - Plots


# Pippi de Bree
# gcd2056@nyu.edu
# 2026-05-14


library(tidyverse)



# Function to summarize convergence 


sims_conv <- readRDS("sims_output/point5_sims/hard_seqs_train_conv_tracker.rds")




  
get_conv <- function(conv_list, link, method, nSoftClusts) {

  all_vecs <- lapply(sims_conv[[link]], function(set_data) {
      set_data[[method]][,nSoftClusts]
    })
  
  combined_vec <- do.call(rbind, all_vecs)
  return(mean(combined_vec))
}

make_full_conv_summary <- function(sim_conv_data) {

  link_types <- c("concord", "indep_semi", "indep", "semi_concord", "semi") 
  methods <- c("om_trate", "om_slog", "lcs")
  nSoft <- 4

  conv_frame <- as.data.frame(expand.grid(link_types, methods, 2:nSoft))
  names(conv_frame) <- c("link", "method", "clusts")

  conv_frame <- conv_frame |> rowwise() |> mutate(prop_converged = get_conv(sim_conv_data, link, method, clusts))

  return(conv_frame)

}



hard_sims_conv <- readRDS("sims_output/point5_sims/hard_seqs_train_conv_tracker.rds")
very_hard_sims_conv <- readRDS("sims_output/point5_sims/very_hard_seqs_train_conv_tracker.rds")
med_sims_conv <- readRDS("sims_output/point5_sims/med_seqs_train_conv_tracker.rds")
dar_med_sims_conv <- readRDS("sims_output/point5_sims/dar_med_seqs_train_conv_tracker.rds")
easy_sims_conv <- readRDS("sims_output/point5_sims/easy_seqs_train_conv_tracker.rds")



hard_train_conv <- make_full_conv_summary(hard_sims_conv)
very_hard_train_conv <- make_full_conv_summary(very_hard_sims_conv)
med_train_conv <- make_full_conv_summary(med_sims_conv)
dar_med_train_conv <- make_full_conv_summary(dar_med_sims_conv)
easy_train_conv <- make_full_conv_summary(easy_sims_conv)





# Have to limit the data such that we make sure there is no issue with the 
# summarizing 
# should also report somewhere how many of each of the soft clustering are actually calculated in the average 

# Transform the list and map the names
plot_data <- map_df(names(mse.sims), function(name) {
  matrix_data <- mse.sims[[name]]
  
  # Calculate row means - only works with n > 1
  row_means <- sqrt(rowMeans(matrix_data))
  
  # Create the data frame
  tibble(
    Method = factor(method_names, levels = method_names), # Keeps them in your specific order
    Mean_MSE = row_means,
    Simulation_Type = name
  )
})



# Create the plot
plot_name <- paste0(sim_type, "plot.pdf")
  
pdf(plot_name,width=9,height=7)
ggplot(plot_data, aes(x = Method, y = Mean_MSE, color = Simulation_Type, group = Simulation_Type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) + # Makes text a bit more readable
  labs(
    title = "Mean RMSE per Method - Hard Difficulty",
    x = "Method",
    y = "Average RMSE (of Best Number of Components)",
    color = "Simulation Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Tilts labels if names are long
dev.off()
  
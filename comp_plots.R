# Code for MVAD Sequence Work Plots

# Pippi de Bree
# gcd2056@nyu.edu

library(tidyverse)


# 1. Load the data
# This is a list of matrices where rows = folds, cols = components/clusters
cov_list <- readRDS("linear_output/linear_job_coverage/LinearFullCoverages.rds")

# 2. Process the list into a long-format data frame
# We calculate the mean across the 5 rows (folds) for each column
cov_long <- map_df(names(cov_list), function(method_name) {
  mat <- cov_list[[method_name]]
  
  # Calculate mean coverage across the 5 folds (rows)
  # We use colMeans and then create a data frame
  means <- colMeans(mat, na.rm = TRUE)
  
  data.frame(
    Method = method_name,
    Index = 1:length(means),
    Coverage = means
  )
})


# 1. Process the data into the long format
cov_long <- map_df(names(cov_list), function(method_name) {
  mat <- cov_list[[method_name]]
  means <- colMeans(mat, na.rm = TRUE)
  
  data.frame(
    Method = method_name,
    Index = 1:length(means),
    Coverage = means
  )
}) %>%
  mutate(Method = case_when(
    Method == "om_trate_hard" ~ "OM T-Rate (Hard)",
    Method == "om_trate_soft" ~ "OM T-Rate (Soft)",
    Method == "om_slog_hard"  ~ "OM Slog (Hard)",
    Method == "om_slog_soft"  ~ "OM Slog (Soft)",
    Method == "lcs_hard"      ~ "LCS (Hard)",
    Method == "lcs_soft"      ~ "LCS (Soft)",
    Method == "windows"       ~ "Windows",
    Method == "harm"          ~ "CFDA",
    TRUE ~ Method
  )) %>%
  filter(Index < 16) # Matching your RMSE plotting logic

# 2. Save to PDF
pdf("LinearCoveragePlot.pdf", width = 8, height = 6)

ggplot(data = cov_long, 
       aes(x = Index, y = Coverage, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  # Often useful to see where 95% coverage lies
  geom_hline(yintercept = 0.95, linetype = "dashed", alpha = 0.5) +
  labs(title = "Linear Regression Coverage by Method",
       subtitle = "Averaged across 5 Folds",
       x = "Components/Clusters",
       y = "Coverage Rate",
       color = "Method") +
  theme_minimal()

dev.off()


# 1. Load the MPIW data
# Structure: List of matrices (Rows = Folds, Cols = Components)
mpiw_list <- readRDS("linear_output/linear_job_coverage/LinearFullMPIW.rds")

# 2. Process the data into a long-format data frame
mpiw_long <- map_df(names(mpiw_list), function(method_name) {
  mat <- mpiw_list[[method_name]]
  
  # Calculate mean interval width across the 5 folds
  means <- colMeans(mat, na.rm = TRUE)
  
  data.frame(
    Method = method_name,
    Index = 1:length(means),
    MPIW = means
  )
}) %>%
  mutate(Method = case_when(
    Method == "om_trate_hard" ~ "OM T-Rate (Hard)",
    Method == "om_trate_soft" ~ "OM T-Rate (Soft)",
    Method == "om_slog_hard"  ~ "OM Slog (Hard)",
    Method == "om_slog_soft"  ~ "OM Slog (Soft)",
    Method == "lcs_hard"      ~ "LCS (Hard)",
    Method == "lcs_soft"      ~ "LCS (Soft)",
    Method == "windows"       ~ "Windows",
    Method == "harm"          ~ "CFDA",
    TRUE ~ Method
  )) %>%
  # Filter to Index < 16 to stay consistent with your RMSE and Coverage plots
  filter(Index < 16)

# 3. Save to PDF
pdf("LinearMPIWPlot.pdf", width = 8, height = 6)

ggplot(data = mpiw_long, 
       aes(x = Index, y = MPIW, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  labs(title = "Linear Regression Interval Width (MPIW) by Method",
       subtitle = "Averaged across 5 Folds (Lower is usually better)",
       x = "Components/Clusters",
       y = "Mean Prediction Interval Width",
       color = "Method") +
  theme_minimal()

dev.off()





cv_lin <- readRDS("cv_outputs/CV_LinearComp.rds")
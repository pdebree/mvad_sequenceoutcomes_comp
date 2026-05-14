# Code for MVAD Sequence Work Plots

# Pippi de Bree
# gcd2056@nyu.edu

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
source("seqout_utils.R")



# Linear Competitition Plots - For CVs 

linear_comp <- readRDS("cv_outputs/CV_LinearComp.rds")
nonlinear_comp <- readRDS("cv_outputs/CV_NonLinearComp.rds")


# MSE Plot
get_method_metrics <- function(data_list, method_name, method_labs) {

  # pull out cross validated mses for given method
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)

  method_cov <- lapply(data_list$cov, function(x) x[[method_name]])
  combined_cov <- simplify2array(method_cov)
  
  method_mpiw <- lapply(data_list$mpiw, function(x) x[[method_name]])
  combined_mpiw <- simplify2array(method_mpiw)

  # find means across the folds (of all cross )
  mean_mse <- apply(combined_mse, 2, mean, na.rm = TRUE)
  mean_cov <- apply(combined_cov, 2, mean, na.rm = TRUE)
  mean_mpiw <- apply(combined_mpiw, 2, mean, na.rm = TRUE)

  df <- data.frame(
      Method = method_name, 
      Index = 1:length(mean_mse),
      RMSE   = sqrt(mean_mse), 
      COV    = mean_cov, 
      MPIW   = mean_mpiw,
      stringsAsFactors = FALSE)


  df <- df |> mutate(Method = recode(Method, !!!method_labs))

  return(df)
}


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

lin_metrics <- do.call(rbind, lapply(methods, function(m) get_method_metrics(linear_comp, m, method_labels)))


lin_rmse_long <- lin_metrics |> dplyr::select(Method, Index, RMSE) |> filter(Index < 16)

# Check for best performance by method 
as.data.frame(lin_rmse_long  |> group_by(Method) |> filter(RMSE == safe_min(RMSE)) |> arrange(RMSE))

lin_rmse_long <- lin_rmse_long |> filter(Method %in% c("Windows", "CFDA", "LCS (Soft)", "OM INDELSLOG (Hard)"))

pdf("plots/LinearCompPlot.pdf",width=8,height=6)
ggplot(data = lin_rmse_long, 
       aes(x = Index, y = RMSE, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  scale_x_continuous(breaks = 2:25) +
  labs(title = "Linear Regression Performance",
       x = "Number of Clusters / Components",
       y = "RMSE (Squared Average MSE Across CVs)",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85), 
    legend.background = element_rect(fill = "white", color = "grey80"))
dev.off()




# ***************************
# Non-Linear 


# function to get best mse (considering mtry) for every method/component combination
get_nonlinear_rmse <- function(data_list, method_name) {
  
  # Extract the list of 20 CV 3d arrays (fold, comps, mtry)
  method_list <- lapply(data_list$mse, function(x) x[[method_name]])

  # make into 4d array 
  combined_4d <- simplify2array(method_list)

  # make matrix of means (i.e. number of components and mtry - means across folds*cvs)
  mean_mse_matrix <- apply(combined_4d, c(2, 3), mean, na.rm = TRUE)

  # make vector of best performance and vector of best associated mtry values
  avg_rmse_by_mtry <- sqrt(apply(mean_mse_matrix, 1, safe_min))
  best_mtry <- apply(mean_mse_matrix, 1, safe_which_min)
  
  return(list(rmses = avg_rmse_by_mtry, mtrys = best_mtry))
}

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

nonlin_metric_data <- map_dfr(names(method_labels), function(m) {
  res <- get_nonlinear_rmse(nonlinear_comp, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmse,
    method = method_labels[[m]], 
    mtry = res$mtrys
  )
})


get_cov_comp_mtry <- function(data_list, method, mtry, comps) {

  method_list <- lapply(data_list$cov, function(x) x[[method_name]])
  
  if (is.na(mtry)) {
    folds_covs <- NA
  } else {
    folds_covs <- mean(simplify2array(method_list)[, comps, mtry, ])
  }
  
  return(folds_covs)

}


get_mpiw_comp_mtry <- function(data_list, method, mtry, comps) {
  # gets the average cov for a number of components, with a given mtry, for a specific method 

  method_list <- lapply(data_list$mpiw, function(x) x[[method_name]])
  
  if (is.na(mtry)) {
    folds_mpiw <- NA
  } else {
    folds_mpiw <- mean(simplify2array(method_list)[, comps, mtry, ])
  }

  return(folds_mpiw)
}

nonlin_metric_data <- nonlin_metric_data |> rowwise() |> mutate(
  cov = get_cov_comp_mtry(data_list, method, mtry, comps),
  mpiw = get_mpiw_comp_mtry(data_list, method, mtry, comps)
)



# Look at best performances 
nonlin_metric_data |> group_by(method) |> filter(rmse == safe_min(rmse)) |> arrange(rmse)

# filter to the CFDA, Windows and the Best Performing hard and Soft Clustering 
nonlin_metric_data <- nonlin_metric_data  |> filter(method %in% c("Counts", "CFDA", ))


pdf("plots/NonLinearCompPlot.pdf", width = 8, height = 6)
ggplot(nonlin_rmse_plot_data, aes(x = comps, y = rmse, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  scale_x_continuous(breaks = 2:25) +
  theme_minimal() +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85), 
    legend.background = element_rect(fill = "white", color = "grey80")) + 
  labs(title = "Random Forest Performance",
       x = "Number of Components / Clusters",
       y = "RMSE (Squared Average MSE Across CVs)",
       color = "Method") 
dev.off()





# Make Tables of Convergences 





# Make Convergence Plots 
cov_list <- readRDS("linear_output/linear_job_coverage/LinearFullCoverages.rds")

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






# Variable Importance Plots 
# Using the mtry associated with best performance


best_mtry <- nonlin_rmse_plot_data |> group_by(method) |> filter(rmse == safe_min(rmse)) 

### Data Preparation 
# Data load in and labeling 
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
                "joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")

# seqdef - creates a state sequence object (formulates sequences into labels)
mvad.seq <- seqdef(mvad, 15:50, states=mvad.scodes, labels=mvad.labels)

# Calculate number of months employed in last year. 
mvad_last_year <- mvad[,75:86]
num_month_em_last_year <- apply(mvad_last_year, 1, function(x) length(which(x=="employment")))
mvad_covars <- mvad[3:14] %>% dplyr::select(-Western) #reference group

# Create distance matrices 
dists <- create_dists(data.seq=mvad.seq)


# Var Importance Plots for Clustering - Come back to check which distance does best overall (for hard and soft)
# Do LCS Soft and OM-Trate Hard 


best_no_hard_cl <- best_mtry |> filter(method == "") |> pull(comps)

# LCS (Hard) - dists[[2]] is hard coded as LCS
clusterward_hard <- agnes(dists[[]], diss=TRUE, method="ward")

mvad_hard_cl <- mvad_covars %>% 
  mutate(cluster=factor(cutree(clusterward_hard,k=best_no_hard_cl)), y=num_month_em_last_year) 

# Create dummy columns for all expected clusters - have to do this explicitly to make 
# sure no columns are dropped if empty in the test set
# do not do cluster 1 because we need to not have a linear combination
for(col in paste0("Cluster_", 2:best_no_hard_cl)) {
  mvad_hard_cl[[col]] <- as.integer(mvad_hard_cl$cluster == sub("Cluster_", "", col))
}
mvad_hard_cl <- mvad_hard_cl %>% dplyr::select(-cluster)


# Soft Clustering - OM T-Rate Soft is the best
best_no_soft_cl <- best_mtry |> filter(method == "")  |> pull(comps)

# OM T-Rate (Soft) - Hard Coded in dists[[1]]
clustering_soft <- fanny(dists[[]], k=best_no_soft_cl, memb.exp=1.5, diss=TRUE, maxit = 1000)$membership
colnames(clustering_soft) <- paste0("Cluster_",1:best_no_soft_cl)

# remove first group so not linearly dependent. 
mvad_soft_cl <- cbind(mvad_covars, clustering_soft[,-1]) |> mutate(y = num_month_em_last_year)
  


### Windows 
mvad_states <- mvad[,c(1,15:50)]

# column names for the three years 
y1 <- colnames(mvad_states)[2:13]
y2 <- colnames(mvad_states)[14:25]
y3 <- colnames(mvad_states)[26:37]

# makes year and id combinations with counts for each state
id_year <- mvad_states %>% pivot_longer(cols=colnames(mvad_states[2:37]), names_to="output") %>% 
    mutate(year = case_when(output %in% y1 ~ 1, output %in% y2 ~ 2, TRUE ~ 3)) %>% count(id, year, value)

# make full table with 16 columns (one for each year/id combo) + 1 for id.
mvad_states_wide <- id_year %>% 
  pivot_wider(id_cols=c(id), 
              names_from=c(value,year), values_from=n, values_fill = 0) %>% dplyr::select(-id)

# do PCA and select number of windows that comes from bst 
pca_windows <- prcomp(x=mvad_states_wide, center=TRUE, scale=TRUE)
windows_pcs <- pca_windows$x[,1:best_no_windows]
colnames(windows_pcs) <- paste0("PC_",1:best_no_windows)

best_no_windows <- best_mtry |> filter(method == "Counts")  |> pull(comps)

mvad_windows <- cbind(num_month_em_last_year, mvad_covars, windows_pcs)
colnames(mvad_windows)[1] <- "y"



# CFDA 
M <- 36-1 #number of initial months-1

mvad_wide <- mvad[,c(1,15:(15+M))] # Goes only to June 96 
colnames(mvad_wide) <- c("id", as.character(0:M))

# have time and id combinations (so basically only one state column)
mvad_long <- gather(mvad_wide, key = time, value = state, "0":as.character(M))
mvad_long$time <- as.numeric(mvad_long$time)
mvad_long <- mvad_long %>% mutate(state=as.factor(state))
summary_cfd(mvad_long)

# first element is the range overwhich the functions can be evaluated 
# nbasis - number of basis functions (here this is 6 because we have 6 states)
# norder - the order of their degree (here cubic - why cubic?)
basis <- create.bspline.basis(c(0, M), nbasis = 6, norder = 4)


fmca <- compute_optimal_encoding(mvad_long, basis, nCores = 7,verbose=F)
pcs.cfda <- fmca$pc # we get 34 pcs (36 month of "variables")

best_no_harms <- best_mtry |> filter(method == "CFDA")  |> pull(comps)

harmonics <- as_tibble(pcs.cfda[, 1:best_no_harms, drop = FALSE])
colnames(harmonics) <- paste0("Harmonic_",1:best_no_harms)

mvad_harm <- cbind(mvad_covars, harmonics) %>% 
  mutate(y = num_month_em_last_year)


# function to make variable importance plot
var_imp_plot <- function(mvad_data, best_mtry, title_string, subtitle_string = "") {
  # must be a "y" column in the dataframe (which is equivalent to num_month_emp)
    # fit with hard clusters 
  fit.rf <- ranger(y ~ .,
                        data = mvad_data,
                        num.trees = 1000, 
                        mtry = best_mtry,
                        keep.inbag = TRUE,
                        respect.unordered.factors = TRUE,
                        quantreg = TRUE, importance = "impurity")
  
  var_imps <- data.frame(fit.rf$variable.importance) |> rownames_to_column("covar") |> 
    rename(var_imp = fit.rf.variable.importance)


  p <- ggplot(data = var_imps, aes(x=var_imp, y=reorder(covar, var_imp))) + geom_col(fill="steelblue") + 
    labs(x="Variable Importance", y="Covariate", title=title_string, subtitle=subtitle_string)

  return(p)


}



best_mtry_windows <- best_avg_mtry |> filter(method == "Counts")  |> pull(mtry)

# Make Variable Importance Plots
wind_vi <- var_imp_plot(mvad_data = mvad_windows, 
  best_mtry = best_mtry_windows,
  title_string = "Variable Importance For Counts", 
  subtitle_string = paste0("Mtry = ",  best_mtry_windows, ", Count Principal Components = ", best_no_windows)
)

pdf("plots/VarImpCounts.pdf", width = 8, height = 6)
wind_vi
dev.off()



best_mtry_harms <- best_mtry  |> filter(method == "CFDA")  |> pull(mtry)

# Make Variable Importance Plots
harm_vi <- var_imp_plot(mvad_data = mvad_harm, 
  best_mtry = best_mtry_harms,
  title_string = "Variable Importance For CFDA", 
  subtitle_string = paste0("Mtry = ",  best_mtry_harms, ", Harmonics = ", best_no_harms)
)

pdf("plots/VarImpHarms.pdf", width = 8, height = 6)
harm_vi
dev.off()



best_mtry_soft <- best_avg_mtry |> filter(method == )  |> pull(mtry)

# Make Variable Importance Plots
soft_vi <- var_imp_plot(mvad_data = mvad_soft_cl, 
  best_mtry = best_mtry_soft,
  title_string = "Variable Importance For HHHHHHH", 
  subtitle_string = paste0("Mtry = ",  best_mtry_harms, ", Soft Clusters = ", best_no_soft_cl)
)

pdf("plots/VarImpSoftCl.pdf", width = 8, height = 6)
soft_vi
dev.off()

best_mtry_hard <- best_mtry |> filter(method == )  |> pull(avg_mtry)

# Make Variable Importance Plots
hard_vi <- var_imp_plot(mvad_data = mvad_hard_cl, 
  best_mtry = best_mtry_hard,
  title_string = "Variable Importance For )", 
  subtitle_string = paste0("Mtry = ",  best_mtry_hard, ", Hard Clusters = ", best_no_hard_cl)
)

pdf("plots/VarImpHardCl.pdf", width = 8, height = 6)
hard_vi
dev.off()

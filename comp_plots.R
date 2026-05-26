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
library(factoextra)
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

  if (method_name == "demos") {
    mean_mse <- mean(combined_mse)
    mean_cov <- mean(combined_cov)
    mean_mpiw <- mean(combined_mpiw)

  } else {
    # find means across the folds (of all cross )
    mean_mse <- apply(combined_mse, 2, mean, na.rm = TRUE)
    mean_cov <- apply(combined_cov, 2, mean, na.rm = TRUE)
    mean_mpiw <- apply(combined_mpiw, 2, mean, na.rm = TRUE)
  }

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
             "om_slog_soft", "lcs_hard", "lcs_soft", "windows", "harm", "demos")

method_labels <- c(
  "om_trate_hard" = "OM T-Rate (Hard)", 
  "om_trate_soft" = "OM T-Rate (Soft)", 
  "om_slog_hard"  = "OM INDELSLOG (Hard)", 
  "om_slog_soft"  = "OM INDELSLOG (Soft)",
  "lcs_hard"      = "LCS (Hard)", 
  "lcs_soft"      = "LCS (Soft)", 
  "windows"       = "Counts", 
  "harm"          = "CFDA", 
  "demos" = "Demographics"
)

lin_metrics <- do.call(rbind, lapply(methods, function(m) get_method_metrics(linear_comp, m, method_labels)))

best_lin <- as.data.frame(lin_metrics |> group_by(Method) |> filter(RMSE == safe_min(RMSE)) |> arrange(RMSE))

lin_rmse_long <- lin_metrics |> dplyr::select(Method, Index, RMSE)  |> filter(Index < 20)

# Check for best performance by method 
as.data.frame(lin_rmse_long  |> group_by(Method) |> filter(RMSE == safe_min(RMSE)) |> arrange(RMSE))


lin_rmse_long <- lin_rmse_long |> filter(Method %in% c("Counts", "CFDA", "LCS (Soft)", "OM INDELSLOG (Hard)"))

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
get_nonlinear_metrics <- function(data_list, method_name) {
  
  # Extract the list of 20 CV 3d arrays (fold, comps, mtry)
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)
  method_cov <- lapply(data_list$cov, function(x) x[[method_name]])
  combined_cov <- simplify2array(method_cov)
  method_mpiw <- lapply(data_list$mpiw, function(x) x[[method_name]])
  combined_mpiw <- simplify2array(method_mpiw)
  
  # within every cross validation - find mse average by mtry 
  avg_cv_mse <- apply(combined_mse, c(2, 3, 4), mean, na.rm = TRUE)
  avg_cv_cov <- apply(combined_cov, c(2, 3, 4), mean, na.rm = TRUE)
  avg_cv_mpiw <- apply(combined_mpiw, c(2, 3, 4), mean, na.rm = TRUE)

  # essentially the location maps to the mtry values
  best_mtry_cv <- apply(avg_cv_mse, c(1, 3), safe_which_min)
  modal_mtrys <- apply(best_mtry_cv, 1, safe_mode)

  # the actual min corresponds to this so we will get mses that correspond to the best
  best_avg_mse_comps_cv <- apply(avg_cv_mse, c(1, 3), safe_min)
  # now we can average these mses (picked by mtry) to find mse by comp
  rmse_comps <- sqrt(apply(best_avg_mse_comps_cv, 1, mean, rm.na=TRUE))

  # data frame that has the indexes needed to get the mtrys for each components/cv combination
  best_mtry_indices <- cbind(
    as.vector(row(best_mtry_cv)), 
    as.vector(best_mtry_cv), 
    as.vector(col(best_mtry_cv))
  )

  # copy dimensions of avg_cv_cov 
  dims <- dim(avg_cv_cov)

  best_avg_cov_comps_cv <- matrix(avg_cv_cov[best_mtry_indices], nrow = dims[1], ncol = dims[3])
  best_avg_mpiw_comps_cv <- matrix(avg_cv_mpiw[best_mtry_indices], nrow = dims[1], ncol = dims[3])

  cov_comps <- apply(best_avg_cov_comps_cv, 1, mean)
  mpiw_comps <- apply(best_avg_mpiw_comps_cv, 1, mean)

  return(list(rmses = rmse_comps, cov = cov_comps, mpiw =  mpiw_comps, modal_mtry = modal_mtrys, mtry_dis = best_mtry_cv, method_name = method_name))

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
  res <- get_nonlinear_metrics(nonlinear_comp, m)  
  data.frame(
    comps = seq_along(res$rmses),
    rmse = res$rmse,
    cov = res$cov, 
    mpiw = res$mpiw,
    method = method_labels[[m]], 
    modal_mtry = res$modal_mtry,
    mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis)))
  )
})


rownames(nonlin_metric_data) <- 1:nrow(nonlin_metric_data)

# Look at best performances 
best_nonlin <- as.data.frame(nonlin_metric_data |> group_by(method) |> filter(rmse == safe_min(rmse)) |> arrange(rmse) |> mutate(sd_mtry = sd(unlist(mtry_dist))))

# filter to the CFDA, Windows and the Best Performing hard and Soft Clustering 
nonlin_metric_data_best <- nonlin_metric_data  |> filter(method %in% c("Counts", "CFDA", "OM INDELSLOG (Hard)", "LCS (Soft)"), comps < 20)


pdf("plots/NonLinearCompPlot.pdf", width = 8, height = 6)
ggplot(nonlin_metric_data_best, aes(x = comps, y = rmse, color = method)) +
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


best_mtry <- nonlin_metric_data |> group_by(method) |> filter(rmse == safe_min(rmse), method %in% c("Counts", "CFDA", "OM INDELSLOG (Hard)", "LCS (Soft)")) 

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


best_no_hard_cl <- best_mtry |> filter(method == "OM INDELSLOG (Hard)") |> pull(comps)

# OM SLOG (Hard) - dists[[3]] is hard coded as OM SLOG 
clusterward_hard <- agnes(dists[[3]], diss=TRUE, method="ward")

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
best_no_soft_cl <- best_mtry |> filter(method == "LCS (Soft)")  |> pull(comps)

# OM T-Rate (Soft) - Hard Coded in dists[[1]]
clustering_soft <- fanny(dists[[2]], k=best_no_soft_cl, memb.exp=1.5, diss=TRUE, maxit = 1000)$membership
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


best_no_windows <- best_mtry |> filter(method == "Counts")  |> pull(comps)

# do PCA and select number of windows that comes from bst 
pca_windows <- prcomp(x=mvad_states_wide, center=TRUE, scale=TRUE)
windows_pcs <- pca_windows$x[,1:best_no_windows]
colnames(windows_pcs) <- paste0("PC_",1:best_no_windows)



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


var_imp_plot <- function(mvad_data, best_mtry, title_string, subtitle_string = "", custom_labs = FALSE) {
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


  if (custom_labs == TRUE) {

      cluster_stats <- mvad_data |> 
        summarize(across(starts_with("Cluster"), ~ sum(.x) / n() * 100)) |>
        pivot_longer(everything(), names_to = "covar", values_to = "pct") |> # Use "covar" here
        mutate(pct_label = paste0(round(pct, 2), "%"))

      var_imps <- var_imps |>
        left_join(cluster_stats |> dplyr::select(covar, pct_label), by = "covar") |>
        mutate(cl_label = ifelse(is.na(pct_label), "", pct_label)) |>
        dplyr::select(-pct_label)

      p <- ggplot(data = var_imps, aes(x=var_imp, y=reorder(covar, var_imp))) + geom_col(fill="steelblue") + 
          geom_text(aes(label = cl_label), size = 3, nudge_x = 1, hjust = 0) + 
        labs(x="Variable Importance", y="Covariate", title=title_string, subtitle=subtitle_string)
    } else {
      p <- ggplot(data = var_imps, aes(x=var_imp, y=reorder(covar, var_imp))) + geom_col(fill="steelblue") + 
        labs(x="Variable Importance", y="Covariate", title=title_string, subtitle=subtitle_string)
    }

  return(p)
}


best_mtry_hard <- best_mtry |> filter(method == "OM INDELSLOG (Hard)" )  |> pull(modal_mtry)

# Make Variable Importance Plots
hard_vi <- var_imp_plot(mvad_data = mvad_hard_cl, 
  best_mtry = best_mtry_hard,
  title_string = "Variable Importance For OM INDELSLOG (Hard)", 
  subtitle_string = paste0("Mtry (Modal) = ",  best_mtry_hard, ", Hard Clusters = ", best_no_hard_cl), custom_labs=TRUE)


pdf("plots/VarImpHardCl.pdf", width = 8, height = 6)
hard_vi
dev.off()

best_mtry_windows <- best_mtry |> filter(method == "Counts")  |> pull(modal_mtry)

# Make Variable Importance Plots
wind_vi <- var_imp_plot(mvad_data = mvad_windows, 
  best_mtry = best_mtry_windows,
  title_string = "Variable Importance For Counts", 
  subtitle_string = paste0("Mtry (Modal) = ",  best_mtry_windows, ", Count Principal Components = ", best_no_windows)
)

pdf("plots/VarImpCounts.pdf", width = 8, height = 6)
wind_vi
dev.off()



best_mtry_harms <- best_mtry  |> filter(method == "CFDA")  |> pull(modal_mtry)

# Make Variable Importance Plots
harm_vi <- var_imp_plot(mvad_data = mvad_harm, 
  best_mtry = best_mtry_harms,
  title_string = "Variable Importance For CFDA", 
  subtitle_string = paste0("Mtry (Modal) = ",  best_mtry_harms, ", Harmonics = ", best_no_harms)
)

pdf("plots/VarImpHarms.pdf", width = 8, height = 6)
harm_vi
dev.off()



best_mtry_soft <- best_mtry |> filter(method == "LCS (Soft)")  |> pull(modal_mtry)

# Make Variable Importance Plots
soft_vi <- var_imp_plot(mvad_data = mvad_soft_cl, 
  best_mtry = best_mtry_soft,
  title_string = "Variable Importance For LCS (Soft)", 
  subtitle_string = paste0("Mtry (Modal) = ",  best_mtry_soft, ", Soft Clusters = ", best_no_soft_cl)
)

pdf("plots/VarImpSoftCl.pdf", width = 8, height = 6)
soft_vi
dev.off()





# Sequence Correlation Plots 


# Sequence Performance plots - demos and non demos

seq_lin_nd <- readRDS("cv_outputs/SeqMetsRes/CV_SeqMetsLinear_NoDemos.rds")
seq_lin_d <- readRDS("cv_outputs/SeqMetsRes/CV_SeqMetsLinear.rds")
seq_nonlin_nd <- readRDS("cv_outputs/SeqMetsRes/CV_SeqMetsNonLinear_NoDemos.rds")
seq_nonlin_d <- readRDS("cv_outputs/SeqMetsRes/CV_SeqMetsNonLinear.rds")


# linear process
lin_seq_best <- function(seq_lin_output) {

  combined_mse <- simplify2array(seq_lin_output$seq_mses)
  seq_3d <- simplify2array(combined_mse)

  comp_means <- sqrt(apply(seq_3d, 2, mean))

  best_n_seq_mets_pc <- safe_which_min(comp_means)

  df <- data.frame(
      Index = 1:length(comp_means),
      RMSE   = comp_means, 
      stringsAsFactors = FALSE)
  return(df)
}


nonlin_seq_best <-  function(seq_nonlin_output) {

  
  seqs_list <- simplify2array(seq_nonlin_output$seq_mses)

  # average across folds
  avg_cv_mse <- apply(seqs_list, c(2, 3, 4), mean, na.rm = TRUE)
  
  # essentially the location maps to the mtry values
  best_mtry_cv <- apply(seqs_list, c(1, 2, 4), safe_which_min)
  # find the mode across the 5 folds 
  modal_mtry_cv <- apply(best_mtry_cv, c(2,3), safe_mode)

  # pull comps based on that modal best 
  modal_mtry_indices <- cbind(
    as.vector(row(modal_mtry_cv)), 
    as.vector(modal_mtry_cv), 
    as.vector(col(modal_mtry_cv))
  )

  # copy dimensions of avg_cv_cov 
  dims <- dim(avg_cv_mse)

  best_avg_mse_comps_cv <- matrix(avg_cv_mse[modal_mtry_indices], nrow = dims[1], ncol = dims[3])

  avg_rmse_comps <- sqrt(apply(best_avg_mse_comps_cv, 1, mean))

  df <- data.frame(
      Index = 1:length(avg_rmse_comps),
      RMSE   = avg_rmse_comps, 
      stringsAsFactors = FALSE)
}


seq_lin_d_rmses <- lin_seq_best(seq_lin_d)
seq_lin_nd_rmses <- lin_seq_best(seq_lin_nd)

seq_nonlin_d_rmses <- nonlin_seq_best(seq_nonlin_d)
seq_nonlin_nd_rmses <- nonlin_seq_best(seq_nonlin_nd)

best_nseq_nonlin_d_rmses <- which.min(seq_nonlin_d_rmses$RMSE)
best_nseq_nonlin_nd_rmses <- which.min(seq_nonlin_nd_rmses$RMSE)
best_seq_nonlin_d_rmse <- min(seq_nonlin_d_rmses$RMSE)
best_seq_nonlin_nd_rmse <- min(seq_nonlin_nd_rmses$RMSE)


print("Nonlin SeqMets - with Demographics")
print(best_nseq_nonlin_d_rmses)
print(best_seq_nonlin_d_rmse)

print("Nonlin SeqMets - w/o Demographics")
print(best_nseq_nonlin_nd_rmses)
print(best_seq_nonlin_nd_rmse)

print("Lin SeqMets - with Demographics")
print(best_nseq_lin_d_rmses)
print(best_seq_lin_d_rmse)

print("Lin SeqMets - w/o Demographics")
print(best_nseq_lin_nd_rmses)
print(best_seq_lin_nd_rmse)




# MSE Plot
get_method_metrics_seq_clust <- function(data_list, method_name) {

  # pull out cross validated mses for given method
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)

  method_cov <- lapply(data_list$cov, function(x) x[[method_name]])
  combined_cov <- simplify2array(method_cov)
  
  method_mpiw <- lapply(data_list$mpiw, function(x) x[[method_name]])
  combined_mpiw <- simplify2array(method_mpiw)

  if (method_name == "demos") {
    mean_mse <- mean(combined_mse)
    mean_cov <- mean(combined_cov)
    mean_mpiw <- mean(combined_mpiw)

  } else {
    # find means across the folds (of all cross )
    mean_mse <- apply(combined_mse, 2, mean, na.rm = TRUE)
    mean_cov <- apply(combined_cov, 2, mean, na.rm = TRUE)
    mean_mpiw <- apply(combined_mpiw, 2, mean, na.rm = TRUE)
  }

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


seq_clust_lin_output <- seq_lin_d

get_method_metrics(seq_lin_d$mse)
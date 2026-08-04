# MVAD Competition Plots 

# Pippi de Bree
# gcd2056 [AT] nyu.edu
# 2026-08-02

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



# read in competition results 
# linear_comp <- readRDS("cv_outputs/upped_nsoftclust/CV_LinearComp_w3pc_175_10.rds")
# nonlinear_comp <- readRDS("cv_outputs/upped_nsoftclust/CV_NonLinearComp_Factor_w3pc_175_10.rds")



linear_comp <- readRDS("cv_outputs/comp_w_seqmets_ccc/MVADLinearComp_Results.rds")
nonlinear_comp <- readRDS("cv_outputs/comp_w_seqmets_ccc/CV_NonLinearComp_Factor_w3pc_r_checks_more_soft.rds")



methods <- c("om_slog_hard", "om_slog_soft", "windows", "harm", "om_slog_hard_3pc", 
            "om_slog_soft_3pc", "demos", "windows_3pc", "harm_3pc")

method_labels <- c(
  "om_slog_hard"  = "OM INDELSLOG (Hard)", 
  "om_slog_soft"  = "OM INDELSLOG (Soft)",
  "windows"       = "Counts", 
  "harm"          = "CFDA",
  "om_slog_hard_3pc"  = "OM INDELSLOG (Hard) + 3PC", 
  "om_slog_soft_3pc"  = "OM INDELSLOG (Soft) + 3PC",
  "windows_3pc" = "Counts + 3PC", 
  "harm_3pc" = "CFDA + 3PC",
  "demos"="Demographics")


# Create dataframe of linear competition 
lin_metrics <- do.call(rbind, lapply(methods, function(m) get_method_metrics(linear_comp, m, method_labels))) 


# Print best results by method
best_lin <- as.data.frame(lin_metrics |> group_by(Method) |> filter(RMSE == safe_min(RMSE)) |> arrange(RMSE))
print(best_lin)


lin_metrics <- lin_metrics |> 
  filter(Method %in% c("Counts", "CFDA", "OM INDELSLOG (Soft)", "OM INDELSLOG (Hard)")) |> 
  mutate(Method = recode_values(
    Method,
    "OM INDELSLOG (Soft)" ~ "Soft Clust.",
    "OM INDELSLOG (Hard)" ~ "Hard Clust.",
    default = Method
  ))

# Create plot of linear competition 
pdf("plots/LinearCompPlot.pdf",width=8,height=6)
ggplot(data = lin_metrics, 
       aes(x = Index, y = RMSE, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  scale_x_continuous(breaks = 2:25) +
  labs(title = "Linear Regression Performance",
       x = "Number of Clusters / Components",
       y = "RMSE",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85), 
    legend.background = element_rect(fill = "white", color = "grey80"))
dev.off()





# Get non-linear metric performance
nonlin_metrics <- get_nonlinear_metrics(nonlinear_comp, method_labels)

# Look at best performances 
best_nonlin <- as.data.frame(nonlin_metrics |> group_by(method) |> filter( rmse == safe_min(rmse)) |> arrange(rmse) |> mutate(sd_mtry = sd(unlist(mtry_dist))))

# filter to the CFDA, Windows and the Best Performing hard and Soft Clustering 
nonlin_metrics  <- nonlin_metrics  |> filter(method %in% c("Counts", "CFDA", "OM INDELSLOG (Hard)", "OM INDELSLOG (Soft)", "OM INDELSLOG (Hard) + 3PC", "OM INDELSLOG (Soft) + 3PC" )) |> 
  mutate(method = recode_values(
    method,
    "OM INDELSLOG (Soft)" ~ "Soft Clust.",
    "OM INDELSLOG (Hard)" ~ "Hard Clust.",
    default = method
  ))



pdf("plots/NonLinearCompPlot.pdf", width = 8, height = 6)
ggplot(nonlin_metrics |> filter(!method %in% c( "OM INDELSLOG (Soft)", "OM INDELSLOG (Hard) + 3PC", "OM INDELSLOG (Soft) + 3PC")), aes(x = comps, y = rmse, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  scale_x_continuous(breaks = 2:25) +
  theme_minimal() +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.75), 
    legend.background = element_rect(fill = "white", color = "grey80")) + 
  labs(title = "Random Forest Performance",
       x = "Number of Components / Clusters",
       y = "RMSE",
       color = "Method") 
dev.off()



### Data Preparation for variance explained plots
data(mvad)
mvad_last_year <- mvad[,75:86]
num_month_em_last_year <- apply(mvad_last_year, 1, function(x) length(which(x=="employment")))
mvad_covars <- mvad[3:14] %>% dplyr::select(-Western) #reference group

# Create distance matrices 
mvad.seq <- create_seq_data(mvad)
dists <- create_dists(mvad.seq)

# Create hierarchical clustering for the whole dataset.
clusterward_slog <- agnes(dists[[3]], diss=TRUE, method="ward")
silhou_plot <- data.frame(nClusts=2:25, sil_width=2:25, ch_index=2:25)

for (i in 2:25) {
  silhou_plot[[i-1, "sil_width"]] <- silhouette(cutree(clusterward_slog, k=i), dists[[1]]) |> 
    as.data.frame() |> pull(sil_width) |> mean()
  silhou_plot[[i-1, "ch_index"]] <- calinhara(dists[[1]], cutree(clusterward_slog, k=i))
}


# put silhouette width and cd index on the same plot
silhou_plot$ch_index <- silhou_plot$ch_index / 1000

silhou_plot_long <- silhou_plot |> pivot_longer(cols = c(ch_index, sil_width), names_to="measurement", values_to="value")

pdf("plots/HardClustGOF.pdf",width=8,height=6)
ggplot(data = silhou_plot_long, aes(x = nClusts, y = value, color = measurement, group = measurement)) + 
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2:25) +
  scale_color_manual(
    values = c("sil_width" = "steelblue", "ch_index" = "darkorange"),
    labels = c("sil_width" = "Average Silhouette Width", "ch_index" = "Calinski–Harabasz Index / 1000")
  ) +
  labs(
    x = "Number of Clusters", 
    y = "Metric Value", 
    title = "Hard Clusters - Goodness of Fit",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    legend.position = "inside", 
    legend.position.inside = c(0.80, 0.1),
    legend.background = element_rect(fill = "white", color = "grey80")
  )
dev.off()



# Soft Clustering Variance Explained 
soft_silhou_plot <- data.frame(nClusts=2:25, sil_width=2:25, ch_index=2:25)

for (i in 2:25) {
  clustering_soft <- fanny(dists[[3]], 
                        k=i, memb.exp=1.75, diss=TRUE, maxit = 2000, tol=1e-10)
  
  soft_silhou_plot[[i-1, "sil_width"]] <- mean(clustering_soft$silinfo$clus.avg.widths)
  soft_silhou_plot[[i-1, "ch_index"]] <- calinhara(dists[[1]], as.vector(apply(clustering_soft$membership, 1, which.max)))
}


# put silhouette width and cd index on the same plot
soft_silhou_plot$ch_index <- soft_silhou_plot$ch_index / 1000

soft_silhou_plot_long <- soft_silhou_plot |> pivot_longer(cols = c(ch_index, sil_width), names_to="measurement", values_to="value")

pdf("plots/SoftlustGOF.pdf",width=8,height=6)
ggplot(data = soft_silhou_plot_long, aes(x = nClusts, y = value, color = measurement, group = measurement)) + 
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2:25) +
  scale_color_manual(
    values = c("sil_width" = "steelblue", "ch_index" = "darkorange"),
    labels = c("sil_width" = "Average Silhouette Width", "ch_index" = "Calinski–Harabasz Index / 1000")
  ) +
  labs(
    x = "Number of Clusters", 
    y = "Metric Value", 
    title = "Soft Clusters - Goodness of Fit",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    legend.position = "inside", 
    legend.position.inside = c(0.80, 0.8),
    legend.background = element_rect(fill = "white", color = "grey80")
  )
dev.off()




# Windows Components
mvad_states <- mvad[,c(1,15:50)]
mvad_states_wide <- create_year_state_counts(mvad)
mvad_windows <- cbind(mvad_covars, mvad_states_wide)

pca_windows_train <- prcomp(x=mvad_windows[12:27], center=TRUE, scale=TRUE)

pdf("plots/WindowsEigenPlot.pdf",width=8,height=6)
fviz_eig(pca_windows_train, choice = "eigenvalue", ncp = 16, , geom = "line") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  labs(x="Number of Principal Components", title="Scree Plot for Counts") + theme_minimal()

dev.off()


# PCA loadings 

year_state_labels <- list(
"school_1" = "School - Year 1",   
"school_2" = "School - Year 2",    
"school_3" = "School - Year 3",  
"training_1" = "Training - Year 1",
"training_2" = "Training - Year 2",
"training_3" = "Training - Year 3",
"FE_1" = "Further Education - Year 1",     
"FE_2" = "Further Education - Year 2",           
"FE_3" = "Further Education - Year 3",     
"HE_3" = "Higher Education - Year 3",         
"employment_1" = "Employment - Year 1",
"employment_2" = "Employment - Year 2",
"employment_3" = "Employment - Year 3",
"joblessness_1" = "Joblessness - Year 1", 
"joblessness_2" = "Joblessness - Year 2", 
"joblessness_3" = "Joblessness - Year 3"
)


wind_pca_loadings <- as.data.frame(pca_windows_train$rotation[, 1:3]) |> rownames_to_column("metric") |> 
  pivot_longer(cols = 2:4, names_to="Component", values_to="loading") |>  
  mutate(year_state = recode(metric, !!!year_state_labels), year_state = factor(year_state, levels= rev(year_state_labels)))


# split the year_state variable into two so we can facet by loading with year as the x axis 
wind_pca_loadings <- separate_wider_delim(wind_pca_loadings, 
                                 cols = year_state, 
                                 delim = " - ", 
                                 names = c("State", "Year"))
wind_pca_loadings$Year <- as.integer(substring(wind_pca_loadings$Year, 6,6))


wind_pca_loadings <- wind_pca_loadings |> mutate(loading = ifelse(Component == "PC2", loading*-1, loading))

pdf("plots/WindPC_Loadings_byYear.pdf",width=8,height=5)
ggplot(
  data = wind_pca_loadings[order(wind_pca_loadings$State == "Higher Education"), ], 
  aes(x = Year, y = loading, color = State, group = State)
) + 
  geom_line(linewidth = 1.2, na.rm = TRUE) + 
  geom_point() +                            
  geom_segment(
    data = subset(wind_pca_loadings, State == "Higher Education"),
    aes(x = 2.5, xend = 3, y = loading, yend = loading),
    linewidth = 1.2, 
    show.legend = FALSE
  ) +
  labs(
    x = NULL, 
    y = "Loading"
  ) +  
  facet_grid(
    . ~ Component, 
    scales = "free_x",
    labeller = labeller(Component = function(x) paste0("Principal Component ", gsub("[^0-9]", "", x)))
  ) + 
  scale_x_continuous(
    name = "Year",
    breaks = seq(0, 3, by = 1)
  ) + 

  scale_color_brewer(palette = "Accent") +
  theme_grey() + 
  theme(
    legend.position = "right",                    
    legend.title = element_text(size = 12),       
    axis.title = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)), 
    strip.background = element_blank(),           
    strip.text = element_text(hjust = 0, size = 11, face = "plain") 
  )
dev.off()


# CFDA
M <- 36 -1
mvad_long <- create_long_data(mvad)
basis <- create.bspline.basis(c(0, M), nbasis = 6, norder = 4)
fmca <- compute_optimal_encoding(mvad_long, basis, nCores = 7,verbose=F)
cfda_eigenvalues <- data.frame(eigenvalues=fmca$eigenvalues) |> rowid_to_column()

pdf("plots/CFDAEigenPlot.pdf",width=8,height=6)
ggplot(data=cfda_eigenvalues, aes(y=eigenvalues, x=rowid)) + geom_line(width=1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1, 34, 2)) +
  scale_y_continuous(breaks= seq(1,27, 2)) + 
  labs(
    x = "Number of Clusters", 
    y = "Metric Value", 
    title = "OM-INDELSLOG Hard Clusters - Goodness of Fit",
    color = "Metric"
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  theme_minimal() + 
  labs(x="Number of Harmonics", title="Scree Plot for CFDA", y="Eigenvalues") 
dev.off()


# Sequence Metrics PCA 
rmetrics <- create_rmetrics(mvad.seq)
melted_cormat <- round(cor(rmetrics), 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Metric1") %>%
  pivot_longer(-Metric1, names_to = "Metric2", values_to = "value")

met_labels <- list(
  "spells"="Spells", "visited_states"="Visited States", "num_of_trans"="Transitions", 
  "entropy"="Entropy", "dss_subs"="DSS Substitutions", "complexity"="Complexity", "turbulence"="Turbulence", "mean_spell_dur"="Mean Spell Duration", "sd_spell_dur"="Spell Duration Std. Dev",   "badness" = "Badness", "Degrad"="Degradation", "insec" = "Insecurity")


melted_cormat <- melted_cormat |> mutate(Metric1 = recode(Metric1, !!!met_labels), Metric2 = recode(Metric2, !!!met_labels))
melted_cormat$Metric1 <- factor(melted_cormat$Metric1, levels=(met_labels))
melted_cormat$Metric2 <- factor(melted_cormat$Metric2, levels=rev(met_labels))

pdf("plots/SeqMet_Corrs.pdf",width=8,height=6)
ggplot(melted_cormat, aes(x = Metric1, y = Metric2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                      midpoint = 0, limit = c(-1, 1), name="Correlation") +
  geom_vline(xintercept = c(7.5, 9.5), color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = c(3.5, 5.5), color = "black", linetype = "dashed", linewidth = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x='', y='', title='Sequence Metric Correlations')
dev.off()

pca_comps_seq <- prcomp(x=rmetrics, center=TRUE, scale=TRUE)
pdf("plots/SeqMetPCs.pdf",width=8,height=6)
fviz_eig(pca_comps_seq, choice = "eigenvalue", ncp = 16, , geom = "line") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  labs(x="Number of Principal Components", title="Scree Plot for Sequence Metrics") + theme_minimal()
dev.off()


# Just the first three eigenvalues 
seqmets_loading <- as.data.frame(pca_comps_seq$rotation) %>%
  mutate(Variable = rownames(.)) %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "component",
    values_to = "loading"
  ) |> filter(component %in% c("PC1", "PC2", "PC3"))
seqmets_loading$Variable <- factor(seqmets_loading$Variable, labels = met_labels)


pdf("plots/SeqMets_3PCLoadings.pdf",width=6,height=6)
ggplot(seqmets_loading, aes(x = Variable, y = loading, fill = component)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ component, ncol = 1) +  # Forces the 4x1 vertical grid
  labs(
    title = "Sequence Metrics Principal Component Loadings",
    x = "Sequence Metrics",
    y = "Loading"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "gray80", fill = NA),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()



# Variable Importance Plots 
best_mtry <- nonlin_metrics |> group_by(method) |> filter(rmse == safe_min(rmse)) 

# Hard Clustering 
best_no_hard_cl <- best_mtry |> filter(method == "Hard Clust.") |> pull(comps)

# OM SLOG (Hard) - dists[[3]] is hard coded as OM SLOG 
clusterward_hard <- agnes(dists[[3]], diss=TRUE, method="ward")

mvad_hard_cl <- mvad_covars %>% 
  mutate(cluster=factor(cutree(clusterward_hard,k=best_no_hard_cl)), y=num_month_em_last_year) 


# Soft Clustering 
best_no_soft_cl <- best_mtry |> filter(method == "Soft Clust.")  |> pull(comps)

# OM T-Rate (Soft) - Hard Coded in dists[[1]]
clustering_soft <- fanny(dists[[3]], k=best_no_soft_cl, memb.exp=1.75, diss=TRUE, maxit = 2000, tol=1e-10)$membership
colnames(clustering_soft) <- paste0("Cluster_",1:best_no_soft_cl)

mvad_soft_cl <- cbind(mvad_covars, clustering_soft[,-1]) |> mutate(y = num_month_em_last_year)
  
### Windows 
best_no_windows <- best_mtry |> filter(method == "Counts")  |> pull(comps)

# do PCA and select number of windows that comes from bst 
pca_windows <- prcomp(x=mvad_states_wide, center=TRUE, scale=TRUE)
windows_pcs <- pca_windows$x[,1:best_no_windows]
colnames(windows_pcs) <- paste0("PC_",1:best_no_windows)

mvad_windows <- cbind(mvad_covars, windows_pcs) |> mutate(y = num_month_em_last_year)



# CFDA 

best_no_harms <- best_mtry |> filter(method == "CFDA")  |> pull(comps)

harmonics <- as_tibble(fmca$pc[, 1:best_no_harms, drop = FALSE])
colnames(harmonics) <- paste0("Harmonic_",1:best_no_harms)

mvad_harm <- cbind(mvad_covars, harmonics) %>% 
  mutate(y = num_month_em_last_year)




# Hard Clustering 
best_no_hard_cl_3pc <- best_mtry |> filter(method == "OM INDELSLOG (Hard) + 3PC") |> pull(comps)

mvad_hard_cl_3pc <- mvad_covars %>% 
  mutate(cluster=factor(cutree(clusterward_hard,k=best_no_hard_cl_3pc )), y=num_month_em_last_year) 

# Soft Clustering 
best_no_soft_cl_3pc <- best_mtry |> filter(method == "OM INDELSLOG (Soft) + 3PC")  |> pull(comps)

# OM T-Rate (Soft) - Hard Coded in dists[[1]]
clustering_soft_3pc <- fanny(dists[[3]], k=best_no_soft_cl_3pc, memb.exp=1.75, diss=TRUE, maxit = 2000, tol=1e-10)$membership
colnames(clustering_soft_3pc) <- paste0("Cluster_",1:best_no_soft_cl_3pc)

mvad_soft_cl_3pc <- cbind(mvad_covars, clustering_soft[,-1]) |> mutate(y = num_month_em_last_year)

pca_loadings <- pca_comps_seq$x[,1:3]
colnames(pca_loadings) <- paste0("SeqMet", colnames(pca_loadings))

# Add in sequence metrics to the hard and soft data frames 
mvad_hard_cl_3pc <- cbind(mvad_hard_cl_3pc, pca_loadings)
mvad_soft_cl_3pc <- cbind(mvad_soft_cl_3pc, pca_loadings)

# Pull modal mtry values for each method
best_mtry_hard <- best_mtry |> filter(method == "Hard Clust." )  |> pull(modal_mtry)
best_mtry_soft <- best_mtry |> filter(method == "Soft Clust.")  |> pull(modal_mtry)
best_mtry_harms <- best_mtry  |> filter(method == "CFDA")  |> pull(modal_mtry)
best_mtry_windows <- best_mtry |> filter(method == "Counts")  |> pull(modal_mtry)
best_mtry_hard_3pc <- best_mtry |> filter(method == "OM INDELSLOG (Hard) + 3PC" )  |> pull(modal_mtry)
best_mtry_soft_3pc <- best_mtry |> filter(method == "OM INDELSLOG (Soft) + 3PC")  |> pull(modal_mtry)



# Make Variable Importance Plots
hard_vi <- var_imp_plot(mvad_data = mvad_hard_cl, 
  best_mtry = best_mtry_hard,
  title_string = "Hard Clustering Variable Importance", 
  subtitle_string = paste0("No. Clusters = ", best_no_hard_cl, ", Mtry = ",  best_mtry_hard))


pdf("plots/VarImpHardCl.pdf", width = 8, height = 6)
hard_vi
dev.off()


wind_vi <- var_imp_plot(mvad_data = mvad_windows, 
  best_mtry = best_mtry_windows,
  title_string = "Counts Variable Importance", 
  subtitle_string = paste0("No. Principal Components = ", best_no_windows,  ", Mtry = ",  best_mtry_windows )
)

pdf("plots/VarImpCounts.pdf", width = 8, height = 6)
wind_vi
dev.off()


harm_vi <- var_imp_plot(mvad_data = mvad_harm, 
  best_mtry = best_mtry_harms,
  title_string = "CFDA Variable Importance", 
  subtitle_string = paste0("No. Harmonics = ", best_no_harms,  ", Mtry = ",  best_mtry_harms) 
)

pdf("plots/VarImpHarms.pdf", width = 8, height = 6)
harm_vi
dev.off()


soft_vi <- var_imp_plot(mvad_data = mvad_soft_cl, 
  best_mtry = best_mtry_soft,
  title_string = "Soft Clustering Variable Importance", 
  subtitle_string = paste0("No. Clusters = ", best_no_soft_cl, ", Mtry = ",  best_mtry_soft)
)

pdf("plots/VarImpSoftCl.pdf", width = 8, height = 6)
soft_vi
dev.off()


hard_3pc_vi <- var_imp_plot(mvad_data = mvad_hard_cl_3pc, 
  best_mtry = best_mtry_hard_3pc,
  title_string = "Hard Clustering w/ Three Sequence Metric Principal Components Variable Importance", 
  subtitle_string = paste0("No. Clusters = ", best_no_hard_cl, ", Mtry = ",  best_mtry_hard_3pc))


pdf("plots/VarImpHardCl_3PC.pdf", width = 8, height = 6)
hard_3pc_vi
dev.off()


soft_3pc_vi <- var_imp_plot(mvad_data = mvad_soft_cl_3pc, 
  best_mtry = best_mtry_soft_3pc,
  title_string = "Soft Clustering w/ Three Sequence Metric Principal Components Variable Importance", 
  subtitle_string = paste0("No. Clusters = ", best_no_soft_cl, ", Mtry = ",  best_mtry_soft_3pc)
)

pdf("plots/VarImpSoftCl_3PC.pdf", width = 8, height = 6)
soft_3pc_vi
dev.off()



# Aside to look into what that high variable importance soft cluster is 
th <- 0.3

high_soft_cl <- mvad_long |> 
  inner_join(soft_cl_high, by = "id") |> 
  drop_na(state) |> 
  count(cluster, time, state, name = "n") |> 
  # Explicitly call tidyr::complete to avoid TraMineR collision
  tidyr::complete(nesting(cluster, time), state, fill = list(n = 0)) |> 
  mutate(props = n / sum(n), .by = c(cluster, time))

pdf("plots/SoftClusterGroupings.pdf", width = 8, height = 6)
ggplot(high_soft_cl, aes(x = time, y = props, group = state, color = state)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(
    ~ cluster, 
    labeller = as_labeller(function(x) gsub("^cl_", "Cluster ", x))
  ) +
  theme_minimal() +
  labs(
    y = "Proportion",
    x = "Time",
    color = "State",
    title = "State Proportions Over Time by Soft Cluster (With Threshold = 0.3)"
  )
dev.off()





# Create variance importance and interaction plots using the vivid package
hard_vii <- vivi_variable_calcs(mvad_data = mvad_hard_cl, 
  best_mtry = best_mtry_hard)


# Make Variable Importance Plots
soft_vii <- vivi_variable_calcs(mvad_data = mvad_soft_cl, 
  best_mtry = best_mtry_soft)

# Make Variable Importance Plots
wind_vii <- vivi_variable_calcs(mvad_data = mvad_windows, 
  best_mtry = best_mtry_windows)


# Make Variable Importance Plots
harm_vii <- vivi_variable_calcs(mvad_data = mvad_harm, 
  best_mtry = best_mtry_harms
)


hard_3pc_vii <- vivi_variable_calcs(mvad_data = mvad_hard_cl_3pc, 
  best_mtry = best_mtry_hard_3pc)


# Make Variable Importance Plots
soft_3pc_vii <- vivi_variable_calcs(mvad_data = mvad_soft_cl_3pc, 
  best_mtry = best_mtry_soft_3pc)


pdf("plots/VarImpInterHardCl.pdf", width = 8, height = 6)
viviHeatmap(
  mat = hard_vii,
  angle = 45  
) +
  ggtitle("Hard Clustering - Variable Importance and Interaction", 
subtitle = paste0("No. Clusters = ", best_no_hard_cl, ", Mtry = ",  best_mtry_hard)) + 
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 8), 
    axis.text.y = element_text(hjust = 1, size = 8),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt") 
  )
dev.off()





pdf("plots/VarImpIntSoftCl.pdf", width = 8, height = 6)
viviHeatmap(
  mat = soft_vii,
  angle = 45  
) +
  ggtitle("Soft Clustering - Variable Importance and Interaction", 
subtitle = paste0("No. Clusters = ", best_no_soft_cl, ", Mtry = ",  best_mtry_soft)) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 8), 
    axis.text.y = element_text(hjust = 1, size = 8),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt") 
  )
dev.off()



pdf("plots/VarImpIntCounts.pdf", width = 8, height = 6)
viviHeatmap(
  mat = wind_vii,
  angle = 45  
) +
  ggtitle("Counts - Variable Importance and Interaction", subtitle=paste0("No. Principal Components = ", best_no_windows,  ", Mtry = ",  best_mtry_windows)) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 8), 
    axis.text.y = element_text(hjust = 1, size = 8),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt") 
  )
dev.off()


pdf("plots/VarImpIntHarm.pdf", width = 8, height = 6)
viviHeatmap(
  mat = wind_vii,
  angle = 45  
) +
  ggtitle("CFDA - Variable Importance and Interaction", subtitle=paste0("No. Principal Components = ", best_no_harms,  ", Mtry = ",  best_mtry_harms)) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 8), 
    axis.text.y = element_text(hjust = 1, size = 8),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt") 
  )
dev.off()

pdf("plots/VarImpIntHard_3PC.pdf", width = 8, height = 6)
viviHeatmap(
  mat = hard_3pc_vii,
  angle = 45  
) +
  ggtitle("Hard Clustering w/ 3 Sequence Metric Principal Components \nVariable Importance and Interaction", subtitle=paste0("No. Harmonics = ", best_no_hard_cl,  ", Mtry = ",  best_mtry_hard_3pc)) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 8), 
    axis.text.y = element_text(hjust = 1, size = 8),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt") 
  )
dev.off()



pdf("plots/VarImpIntSoft_3PC.pdf", width = 8, height = 6)
viviHeatmap(
  mat = soft_3pc_vii,
  angle = 45  
) +
  ggtitle("Soft Clustering w/ 3 Sequence Metric Principal Components \nVariable Importance and Interaction", subtitle=paste0("No. Harmonics = ", best_no_soft_cl,  ", Mtry = ",  best_mtry_soft_3pc)) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 8), 
    axis.text.y = element_text(hjust = 1, size = 8),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt") 
  )
dev.off()



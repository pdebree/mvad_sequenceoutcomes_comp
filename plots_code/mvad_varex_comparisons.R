# Pippi de Bree
# Code to make visualizations of variance explained by each method. 
# Methods: 

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
library(factoextra)
library(fpc)
library(ggplot2)
library(RColorBrewer)
source("seqout_utils.R")


# Decision to use the whole dataset - check with Marc

# make plots for the best number of components for each method 

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



# pdf("plots/HardClustSilhouette.pdf",width=8,height=6)
# ggplot(data=silhou_plot, aes(x=nClusts, y=sil_width)) + geom_col(fill="steelblue") + 
#   labs(x="Number of Clusters", y="Silhouette Width", title="OM-INDELSLOG Hard Clusters by Average Silhouette Width")
# dev.off()

# pdf("plots/HardClusterCH.pdf",width=8,height=6)
# ggplot(data=silhou_plot, aes(x=nClusts, y=ch_index)) + geom_col(fill="steelblue") + 
#   labs(x="Number of Clusters", y="Calinski-Harabasz Distance", title="OM-INDELSLOG Hard Clusters by Calinski-Harabasz Distance")
# dev.off()

# Windows Components
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

mvad_windows <- cbind(mvad_covars, mvad_states_wide)

pca_windows_train <- prcomp(x=mvad_windows[12:27], center=TRUE, scale=TRUE)

# pdf("plots/WindowsVarEx.pdf",width=8,height=6)
# fviz_eig(pca_windows_train, choice = "variance", ncp = 16, addlabels = TRUE) + 
#     labs(x="Number of Components", title="Scree Plot of Variance Explained for Windows Counts") 
#   theme_minimal() +
# dev.off()

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
  geom_line(linewidth = 1.2, na.rm = TRUE) + # Will now draw Higher Ed on top
  geom_point() +                             # Will now draw Higher Ed on top
  
  # Adds a small horizontal bar only for the Higher Education point
  geom_segment(
    data = subset(wind_pca_loadings, State == "Higher Education"),
    aes(x = 2.5, xend = 3, y = loading, yend = loading),
    linewidth = 1.2, 
    show.legend = FALSE
  ) +

  # 1. Remove the global x label from labs so it doesn't span across the center bottom
  labs(
    x = NULL, 
    y = "Loading"
  ) +  
  
  # 2. Use facet_grid with free_x to isolate the axes completely per plot
  facet_grid(
    . ~ Component, 
    scales = "free_x",
    labeller = labeller(Component = function(x) paste0("Principal Component ", gsub("[^0-9]", "", x)))
  ) + 
  
  # 3. Use 'name' here to repeat the "Year" title at the bottom of EVERY column
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
    axis.title.x = element_text(margin = margin(t = 10)), # Gives spacing to the repeated X titles
    strip.background = element_blank(),           
    strip.text = element_text(hjust = 0, size = 11, face = "plain") 
  )
dev.off()






# CFDA - put eigenvalues directly into this

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


# plotEigenvalues(fmca) + geom_hline(yintercept = 1, color="red") + 
#   scale_fill_manual("steelblue") +
#   scale_x_continuous(n.breaks = 14) + 
#   scale_y_continuous(n.breaks = 35) + 
#   labs(x="Number of Harmonics", title="Eigenvalues Associated with Harmonics")


# pdf("plots/CFDAVarExPlot.pdf",width=8,height=6)
# plotEigenvalues(fmca, normalize=TRUE, cumulative=TRUE) + 
#   scale_fill_manual("steelblue") +
#   scale_x_continuous(n.breaks = 10) + 
#   scale_y_continuous(n.breaks = 10) + 
#   labs(x="Number of Harmonics", title="Variance Explained Associated with Harmonics", 
# y="Cumulative Proportion of Variance Explained")
# dev.off()


# soft 
soft_silhou_plot <- data.frame(nClusts=2:13, sil_width=2:13, ch_index=2:13)

for (i in 2:13) {
  clustering_soft <- fanny(dists[[2]], 
                        k=i, memb.exp=1.5, diss=TRUE, maxit = 1000)
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
    legend.position.inside = c(0.80, 0.1),
    legend.background = element_rect(fill = "white", color = "grey80")
  )
dev.off()




# pdf("plots/SoftClustSilhouette.pdf",width=8,height=6)
# ggplot(data=soft_silhou_plot, aes(x=nClusts, y=sil_width)) + geom_col(fill="steelblue") + 
#   labs(x="Number of Clusters", y="Silhouette Width", title="LCS Soft Clusters by Average Silhouette Width")
# dev.off()

# pdf("plots/SoftClustCH.pdf",width=8,height=6)
# ggplot(data=soft_silhou_plot, aes(x=nClusts, y=ch_index)) + geom_col(fill="steelblue") + 
#   labs(x="Number of Clusters", y="Calinski-Harabasz Distance", title="LCS Soft Clusters by Calinski-Harabasz Distance")
# dev.off()












# Sequence Metrics PCA 


mvad_states <- mvad[,15:50]

# encodings for the statebadness 
st_alphabet <- alphabet(mvad.seq) 
# Example: If alphabet is "A", "B", "C"
st_prec_values <- c(1, -1, -1, 2, -1, -1) # A=1 (low badness), C=3 (high badness)
names(st_prec_values) <- st_alphabet

mvad_rmetrics <- mvad_states %>% mutate(
  spells=seqindic(mvad.seq, "dlgth")$Dlgth, 
  visited_states=seqindic(mvad.seq, "visited")$Visited, 
  num_of_trans = seqindic(mvad.seq,"trans")$Trans, 
  mean_spell_dur = seqindic(mvad.seq,"meand")$MeanD, 
  # pedantic - this one pulled out a "seqivardur" "numeric" datatype (not sure 
  # how it did both, so I have to force it to be numeric)
  sd_spell_dur = as.numeric(seqindic(mvad.seq,"dustd")$Dustd), 
  # Diversity I
  entropy = seqindic(mvad.seq,"entr")$Entr, 
  # more interested in states than spells (see paper)
  dss_subs = seqindic(mvad.seq,"nsubs")$Nsubs, 
  complexity = seqindic(mvad.seq,"cplx")$Cplx, 
  # could look at other turbulence measures 
  turbulence = seqindic(mvad.seq,"turb")$Turb, 
  badness = seqibad(seqdata = mvad.seq,stprec = st_prec_values), 
  degradation = seqidegrad(seqdata = mvad.seq, stprec = st_prec_values), 
  insecurity = seqinsecurity(seqdata = mvad.seq, stprec = st_prec_values))

rmetrics <- mvad_rmetrics[,37:48] 

pca_comps_seq <- prcomp(x=rmetrics, center=TRUE, scale=TRUE)
    
pdf("plots/SeqMetPCs.pdf",width=8,height=6)

fviz_eig(pca_comps_seq, choice = "eigenvalue", ncp = 16, , geom = "line") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  labs(x="Number of Principal Components", title="Scree Plot for Sequence Metrics") + theme_minimal()

dev.off()


metric_labels <- list(
"spells"="Number of Spells", 
"visited_states" ="Numner of Visited States", 
"num_of_trans"="Number of Transitions",  
"mean_spell_dur"="Mean Spell Duration",  
"sd_spell_dur"="Standard Deviation of Spell Duration",
"entropy"="Entropy",        
"dss_subs" = "Number of DSS Subsequences",
"complexity" = "Complexity Index", 
"turbulence" = "Turbulence", 
"badness" = "Badness",
"degradation"  = "Degradation Index", 
"insecurity" = "Insecurity Index"    
)

metric_labels <- list(
"spells"="Number of Spells", 
"visited_states" ="Numner of Visited States", 
"num_of_trans"="Number of Transitions",  
"mean_spell_dur"="Mean Spell Duration",  
"sd_spell_dur"="Standard Deviation of Spell Duration",
"entropy"="Entropy",        
"dss_subs" = "Number of DSS Subsequences",
"complexity" = "Complexity Index", 
"turbulence" = "Turbulence", 
"badness" = "Badness",
"degradation"  = "Degradation Index", 
"insecurity" = "Insecurity Index"    
)



# Just the first three eigenvalues 

pca_loadings <- as.data.frame(pca_comps_seq$rotation[, 1:3]) |> rownames_to_column("metric") |> 
  pivot_longer(cols = 2:4, names_to="Component", values_to="loading") |> 
  mutate(metric = recode(metric, !!!metric_labels), metric = factor(metric, levels= rev(metric_labels))) 
 

pdf("plots/SeqMetPC_Loadings.pdf",width=8,height=6)
ggplot(data=pca_loadings, aes(x=loading, y=metric, color=Component)) + 
  geom_point() + 
  labs(x = "Loading", y="Metric", title="Metric Principal Component Loadings") +  
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed")
dev.off()






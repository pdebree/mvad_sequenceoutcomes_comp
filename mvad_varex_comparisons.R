# Pippi de Bree
# Code to make visualizations of variance explained by each method. 
# Methods: 

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
library(factoextra)
library(fpc)
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
clusterward_trate <- agnes(dists[[1]], diss=TRUE, method="ward")
silhou_plot <- data.frame(nClusts=2:25, sil_width=2:25, ch_index=2:25)

for (i in 2:25) {
  silhou_plot[[i-1, "sil_width"]] <- silhouette(cutree(clusterward_trate, k=i), dists[[1]]) |> 
    as.data.frame() |> pull(sil_width) |> mean()
  silhou_plot[[i-1, "ch_index"]] <- calinhara(dists[[1]], cutree(clusterward_trate, k=i))
}



pdf("plots/HardClustSilhouette.pdf",width=8,height=6)
ggplot(data=silhou_plot, aes(x=nClusts, y=sil_width)) + geom_col(fill="steelblue") + 
  labs(x="Number of Clusters", y="Silhouette Width", title="OM-Trate Hard Clusters by Average Silhouette Width")
dev.off()

pdf("plots/HardClusterCH.pdf",width=8,height=6)
ggplot(data=silhou_plot, aes(x=nClusts, y=ch_index)) + geom_col(fill="steelblue") + 
  labs(x="Number of Clusters", y="Calinski-Harabasz Distance", title="OM-Trate Hard Clusters by Calinski-Harabasz Distance")
dev.off()

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

pdf("plots/WindowsVarEx.pdf",width=8,height=6)
fviz_eig(pca_windows_train, choice = "variance", ncp = 16, addlabels = TRUE) + 
    labs(x="Number of Components", title="Scree Plot of Variance Explained for Windows Counts")
dev.off()

pdf("plots/WindowsEigenPlot.pdf",width=8,height=6)
fviz_eig(pca_windows_train, choice = "eigenvalue", ncp = 16) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  labs(x="Number of Components", title="Scree Plot of Eigenvalues for Windows Counts")
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

pdf("plots/CFDAEigenPlot.pdf",width=8,height=6)
plotEigenvalues(fmca) + geom_hline(yintercept = 1, color="red") + 
  scale_fill_manual("steelblue") +
  scale_x_continuous(n.breaks = 10) + 
  scale_y_continuous(n.breaks = 10) + 
  labs(x="Number of Harmonics", title="Eigenvalues Associated with Harmonics")
dev.off()


# soft 
soft_silhou_plot <- data.frame(nClusts=2:13, sil_width=2:13, ch_index=2:13)

for (i in 2:13) {
  clustering_soft <- fanny(dists[[2]], 
                        k=i, memb.exp=1.5, diss=TRUE, maxit = 1000)
  soft_silhou_plot[[i-1, "sil_width"]] <- mean(clustering_soft$silinfo$clus.avg.widths)
  soft_silhou_plot[[i-1, "ch_index"]] <- calinhara(dists[[1]], as.vector(apply(clustering_soft$membership, 1, which.max)))
}

pdf("plots/SoftClustSilhouette.pdf",width=8,height=6)
ggplot(data=soft_silhou_plot, aes(x=nClusts, y=sil_width)) + geom_col(fill="steelblue") + 
  labs(x="Number of Clusters", y="Silhouette Width", title="LCS Soft Clusters by Average Silhouette Width")
dev.off()

pdf("plots/SoftClustCH.pdf",width=8,height=6)
ggplot(data=soft_silhou_plot, aes(x=nClusts, y=ch_index)) + geom_col(fill="steelblue") + 
  labs(x="Number of Clusters", y="Calinski-Harabasz Distance", title="LCS Soft Clusters by Calinski-Harabasz Distance")
dev.off()

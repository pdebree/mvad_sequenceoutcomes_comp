# CFDA Linear Competition 
# First written for Supervised and Unsupervised Machine Learning Practicum in January 2026


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

source("mvad_seqout_functions.R")


# Data Load in 
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
                 "joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")

# seqdef - creates a state sequence object (formulates sequences into labels)
mvad.seq <- seqdef(mvad, 15:50, states=mvad.scodes, labels=mvad.labels)

# Creates distance matrix 
dists <- create_dists(data.seq=mvad.seq)


# Training set up 
folds <-  5
nrecs <- nrow(mvad[1])

# Set up training and test split so it is the same across the competition 
idx.orig <- rep(1:folds,each=floor(nrecs/folds))
if (nrecs %% folds != 0) idx.orig <- c(idx.orig,1:(nrecs %% folds))
idx <- sample(idx.orig) #shuffle

# Creates a random shuffle of the indices to be used in the the train-test split and folds. 
shuffle_index <- sample(1:nrecs, nrecs)
ids <- sort(unique(mvad$id))

# set up data for predictions 
mvad_last_year <- mvad[,75:86]
num_month_em_last_year <- apply(mvad_last_year, 1, function(x) length(which(x=="employment")))
mvad_covars <- mvad[3:14] %>% dplyr::select(-Western) #reference group

maxMtry <- 11
nWindows <- 16 
nClusts <- 25 
nSoftClusts <- 14 
nHarms <- 25



# Arrays for holding outcomes fits 
rmse.cv.harm_rf <- array(NA,c(folds,nHarms, maxMtry))
rmse.cv.windows_rf <- array(NA,c(folds,nWindows,maxMtry))
rmse.cv.om_hard_rf <- array(NA,c(folds,nClusts,maxMtry))
rmse.cv.om_soft_rf <- array(NA,c(folds,nSoftClusts,maxMtry))
rmse.cv.lcs_hard_rf <- array(NA,c(folds,nClusts,maxMtry))
rmse.cv.lcs_soft_rf <- array(NA,c(folds,nSoftClusts,maxMtry))
rmse.cv.om_seq_hard <- array(NA,c(folds,nClusts,maxMtry))
rmse.cv.om_seq_soft <- array(NA,c(folds,nSoftClusts,maxMtry))

### CFDA Cross Validation 

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

mvad_long0 <- mvad_long
ids <- sort(unique(mvad_long$id))


for (i in 1:folds) {  
  cat("Fold Number ",i,"\n")

  # pull data
  test_idx <- idx == i  
  train_idx <- !test_idx
  mvad_long_train <- mvad_long %>% filter(id %in% ids[train_idx])
  mvad_long_test <- mvad_long %>% filter(id %in% ids[test_idx])
  
  # compute encodings for - get warning messages that at least one states not in support of 
  # one basis function (I assume this is because of HE in the earlier years)
  fmca.train <- compute_optimal_encoding(mvad_long_train, basis, nCores = 7,verbose=F)
  pcs.train <- fmca.train$pc # we get 34 pcs (36 month of "variables")
  nComps <- ncol(pcs.train) 
  
  colnames(pcs.train) <- paste0("PC",1:nComps)
  pcs.test <- predict(fmca.train,newdata=mvad_long_test,method="parallel",nCores=7)
  colnames(pcs.test) <- paste0("PC",1:nComps)

  for (j in 1:nHarms) {
    for (k in 1:maxMtry) {
        # create training set based on number of PCs to include
      train_harm <- mvad_covars[train_idx, ] %>% 
        add_column(as_tibble(pcs.train[, 1:j, drop = FALSE])) %>% 
        mutate(num_month_em_last_year = num_month_em_last_year[train_idx])
      test_harm <- mvad_covars[test_idx, ] %>% 
        add_column(as_tibble(pcs.test[, 1:j, drop = FALSE])) %>% 
        mutate(num_month_em_last_year = num_month_em_last_year[test_idx])
      if (k < ncol(train_harm) - 1) {
        rmse.cv.harm_rf[i, j, k] <- train_rmse_rf(train_harm, test_harm, mtry =
                                                    k)
      }
    }
  }
}



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

mvad_windows <- cbind(mvad_covars, mvad_states_wide)


for (i in 1:folds) {

  cat("Fold Number ",i,"\n")
  test_idx <- idx == i  
  train_idx <- !test_idx

  # pull out training and testing for only the windows
  train_windows <- mvad_windows[12:27][train_idx, ]
  test_windows <- mvad_windows[12:27][test_idx, ]
  
  pca_windows_train <- prcomp(x=train_windows, center=TRUE, scale=TRUE)
  train_scores <- pca_windows_train$x
  test_scores <- predict(pca_windows_train, newdata = test_windows) 
  
  # add principal components to demographic data (no longer need original windows)
  train_mvad_windows <- cbind(num_month_em_last_year[train_idx], mvad_windows[1:11][train_idx,], train_scores) 
  test_mvad_windows <- cbind(num_month_em_last_year[test_idx], mvad_windows[1:11][test_idx,], test_scores)
  
  colnames(train_mvad_windows)[1] <- "num_month_em_last_year"
  colnames(test_mvad_windows)[1] <- "num_month_em_last_year"
  
  for (j in 1:nWindows) {
    for (k in 1:maxMtry) {
      rmse.cv.windows_rf[i,j,k] <- train_rmse_rf(train_data=train_mvad_windows[,1:(12+j)], test_data = test_mvad_windows, mtry=k)
    }
  }
}




## OM-TRATE


# go through folds in repeat
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx

  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
  clusterward_hard <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")

  # hard coded for OM-Trate
  for (j in 2:nClusts) {
    
    # Create hard clusters
    hard_cluster_data <- hard_cluster(
      clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, 
      num_month_em_last_year = num_month_em_last_year, train_idx=train_idx, 
      test_idx=test_idx, dist_matrix = dists[[1]])

    train_om_hard <- hard_cluster_data$train_data
    test_om_hard <- hard_cluster_data$test_data
  
    # Create soft clusters
    soft_cluster_data <- soft_cluster(dist_matrix = dists[[1]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year)

    train_om_soft <- soft_cluster_data$train_data
    test_om_soft <- soft_cluster_data$test_data

    for (k in 1:maxMtry) {
      rmse.cv.om_hard_rf[i,j,k] <- train_rmse_rf(train_data=train_om_hard, test_data=test_om_hard, mtry=k)
      
      if (j <= nSoftClusts) {
        rmse.cv.om_soft_rf[i,j,k] <- train_rmse_rf(train_data=train_om_soft, test_data=test_om_soft, mtry=k)
      }
    }
  }
}

# LCS

# go through folds in repeat
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx

  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
  clusterward_hard <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")

  for (j in 2:nClusts) {
    
    cluster_data <- hard_cluster(clusterward = clusterward_hard, nClusts = j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year, train_idx = train_idx, test_idx = test_idx, dist_matrix = dists[[2]])
    
    train_lcs_hard <- cluster_data$train_data
    test_lcs_hard <- cluster_data$test_data

    soft_cluster_data <- soft_cluster(dist_matrix = dists[[2]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year)

    train_lcs_soft <- soft_cluster_data$train_data
    test_lcs_soft <- soft_cluster_data$test_data
    
    for (k in 1:maxMtry) {

      rmse.cv.lcs_hard_rf[i,j,k] <- train_rmse_rf(train_data=train_lcs_hard, test_data=test_lcs_hard, mtry=k)
      if (j <= nSoftClusts) {
        rmse.cv.lcs_soft_rf[i,j,k] <- train_rmse_rf(train_data=train_lcs_soft, test_data=test_lcs_soft, mtry=k)
      }
    }
  }
}

## Clustering + Sequence Metrics 

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


# Training set up 
rmetrics <- mvad_rmetrics[,37:48] 
nSeqPcs <- ncol(rmetrics) - 1


# go through folds in repeat
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx
  
  train_seqmets <- rmetrics[train_idx, ]
  test_seqmets <- rmetrics[test_idx, ]
  
  pca_comps_train <- prcomp(x=train_seqmets, center=TRUE, scale=TRUE)
  
  train_scores <- pca_comps_train$x
  test_scores <- predict(pca_comps_train, newdata = test_seqmets) 
  
  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
  clusterward_hard <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")

  # keep as matrices, for easier combination with covariate data 
  train_k_scores <- train_scores[, 1:nSeqPcs, drop = FALSE]
  test_k_scores <- test_scores[, 1:nSeqPcs, drop = FALSE]
    
    for (j in 2:nClusts) {
      
      cluster_hard <- hard_cluster(clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year, train_idx=train_idx, test_idx = test_idx, dist_matrix = dists[[1]])
      
      train_om_hard <- cbind(cluster_hard$train_data, train_k_scores)
      test_om_hard <- cbind(cluster_hard$test_data, test_k_scores) 
      
      if (j < 16) {
      cluster_soft <- soft_cluster(
        dist_matrix = dists[[1]],
        train_idx = train_idx,
        test_idx = test_idx,
        nClusts = j,
        covars = mvad_covars,
        num_month_em_last_year = num_month_em_last_year)
      
      # add in scores data 
      train_om_soft <- cbind(cluster_soft$train_data, train_k_scores)
      test_om_soft <- cbind(cluster_soft$test_data, test_k_scores)
      }
      for (k in 1:maxMtry) {
        rmse.cv.om_seq_hard[i,j,k] <- train_rmse_rf(
          train_data=train_om_hard, test_data=test_om_hard, mtry=k)
        if (j <= nSoftClusts) {
          rmse.cv.om_seq_soft[i,j,k] <- train_rmse_rf(
          train_data=train_om_soft, test_data=test_om_soft, mtry=k)
        }
      
    }
  }
}


harm_means_rf <- apply(rmse.cv.harm_rf, c(2,3), mean)
windows_means_rf <- apply(rmse.cv.windows_rf, c(2,3), mean)
om_hard_means_rf <- apply(rmse.cv.om_hard_rf, c(2,3), mean)
om_soft_means_rf <- apply(rmse.cv.om_soft_rf, c(2,3), mean)
lcs_hard_means_rf <- apply(rmse.cv.lcs_hard_rf, c(2,3), mean)
lcs_soft_means_rf <- apply(rmse.cv.lcs_soft_rf, c(2,3), mean)
soft_cluster_seqs <- apply(rmse.cv.om_seq_soft, c(2,3), mean)
hard_cluster_seqs <- apply(rmse.cv.om_seq_hard, c(2,3), mean)


# Long form 
windows_long <- pivot_means(windows_means_rf, "Windows", nWindows)
harms_long <- pivot_means(harm_means_rf, "Harmonics", nHarms)
om_hard_long <- pivot_means(om_hard_means_rf, "Clusters", nClusts)
om_soft_long <- pivot_means(om_soft_means_rf, "Clusters", nSoftClusts)
lcs_hard_long <- pivot_means(lcs_hard_means_rf, "Clusters", nClusts)
lcs_soft_long <- pivot_means(lcs_soft_means_rf, "Clusters", nSoftClusts)
soft_seq_long <- pivot_means(soft_cluster_seqs, "Clusters", nSoftClusts)
hard_seq_long <- pivot_means(hard_cluster_seqs, "Clusters", nClusts)


best_windows <- windows_long %>% group_by(Windows) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_harms <- harms_long %>% group_by(Harmonics) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_om_hard <- om_hard_long %>% group_by(Clusters) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_om_soft <- om_soft_long %>% group_by(Clusters) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_lcs_hard <- lcs_hard_long %>% group_by(Clusters) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_lcs_soft <- lcs_soft_long %>% group_by(Clusters) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_soft_seq <- soft_seq_long %>% group_by(Clusters) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup
best_hard_seq <- hard_seq_long %>% group_by(Clusters) %>% drop_na() %>% slice_min(RMSE,n=1) %>% mutate(RMSE = sprintf("%.5f", RMSE)) %>% ungroup


# Create data frame of all data 
best.rmse <- data.frame(matrix(NA, nrow = nClusts, ncol = 1))
colnames(best.rmse) <- c("seq_soft")
best.rmse$seq_soft <- c(NA, best_soft_seq$RMSE, rep(NA, nClusts-nSoftClusts))
best.rmse$seq_hard <- c(NA, best_hard_seq$RMSE)
best.rmse$cfda <- c(NA, best_harms$RMSE)
best.rmse$windows <- c(best_windows$RMSE, rep(NA, nClusts - nWindows))
best.rmse$om_hard <- c(NA, best_om_hard$RMSE)
best.rmse$om_soft <- c(NA, best_om_soft$RMSE, rep(NA, nClusts-nSoftClusts))
best.rmse$lcs_hard <- c(NA, best_lcs_hard$RMSE)
best.rmse$lcs_soft <- c(NA, best_lcs_soft$RMSE, rep(NA, nClusts-nSoftClusts))


best.rmse$index <- 1:nrow(best.rmse)

best.rmse_long <- best.rmse %>%
  pivot_longer(
    cols = -index,          
    names_to = "Method",    
    values_to = "RMSE"      
  ) %>% rename(Component=index)


best.rmse_long$RMSE <- as.numeric(as.character(best.rmse_long$RMSE))

pdf("NonLinearCompPlot.pdf",width=9,height=7)
ggplot(data = best.rmse_long, 
       aes(x = Component, y = RMSE, color = Method,group=Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  labs(title = "Random Forest Tuned Performance by Method",
       x = "Components/Clusters",
       y = "RMSE Value",
       color = "Method") +
  scale_color_discrete(labels = c(cfda = "CFDA", lcs_hard = "Hard Clusters (LCS)", lcs_soft = "Hard Clusters (LCS)", om_hard="Hard Clusters (OM-TRate)", om_soft="Soft Clusters (OM-Trate)",  seq_hard="Hard Clusters + 9 Sequence Metric PCs", seq_soft="Soft Clusters + 9 Sequence Metric PCs", windows="Windows")) + theme_minimal()
dev.off()

write.csv(best.rmse, "NonLinearCompetitionRMSE.csv", row.names=FALSE)







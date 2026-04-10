# CFDA Nonlinear Competition 
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
set.seed(111)


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


nWindows <- 16 
nClusts <- 25
nSoftClusts <- 13
nHarms <- 25
nCovars <- ncol(mvad_covars)
nSeqPcs <- 12
fuzz_soft <- 1.5


# Arrays for holding outcomes fits 
# soft clusters - 2 (because we never look at index=1 clusters and indexing is always from 1)
mse.cv.harm_rf <- array(NA,c(folds,nHarms, nCovars + nHarms - 1 ))
mse.cv.windows_rf <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))
mse.cv.om_trate_hard_rf <- array(NA,c(folds,nClusts, nCovars))
mse.cv.om_trate_soft_rf <- array(NA,c(folds,nSoftClusts,nCovars + nSoftClusts - 2))
mse.cv.om_slog_hard_rf <- array(NA,c(folds,nClusts, nCovars))
mse.cv.om_slog_soft_rf <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
mse.cv.lcs_hard_rf <- array(NA,c(folds,nClusts, nCovars))
mse.cv.lcs_soft_rf <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
mse.cv.rmets_rf <- array(NA,c(folds, nSeqPcs, nCovars + nSeqPcs - 1))


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
    for (k in 1:(floor(sqrt(nCovars + j - 1)))) {
        # create training set based on number of PCs to include
      train_harm <- mvad_covars[train_idx, ] %>% 
        add_column(as_tibble(pcs.train[, 1:j, drop = FALSE])) %>% 
        mutate(num_month_em_last_year = num_month_em_last_year[train_idx])
      test_harm <- mvad_covars[test_idx, ] %>% 
        add_column(as_tibble(pcs.test[, 1:j, drop = FALSE])) %>% 
        mutate(num_month_em_last_year = num_month_em_last_year[test_idx])
      if (k < ncol(train_harm) - 1) {
        mse.cv.harm_rf[i, j, k] <- train_mse_rf(train_harm, test_harm, mtry =
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
    for (k in 1:(floor(sqrt(nCovars + j - 1)))) {
      mse.cv.windows_rf[i,j,k] <- train_mse_rf(
        train_data=train_mvad_windows[,1:(12+j)], test_data = test_mvad_windows, mtry=k)
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
    hard_cluster_data <- hard_cluster_onehot(
      clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, 
      num_month_em_last_year = num_month_em_last_year, train_idx=train_idx, 
      test_idx=test_idx, dist_matrix = dists[[1]])

    train_om_hard <- hard_cluster_data$train_data
    test_om_hard <- hard_cluster_data$test_data
  
    # Create soft clusters
    soft_cluster_data <- soft_cluster(dist_matrix = dists[[1]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year)

    train_om_soft <- soft_cluster_data$train_data
    test_om_soft <- soft_cluster_data$test_data

    # have to minus 2 because of indexing j from 2 (instead of 1 - can't look at 1 soft cluster)
    for (k in 1:(floor(sqrt(nCovars + j - 2)))) {

      mse.cv.om_trate_hard_rf[i,j,k] <- train_mse_rf(train_data=train_om_hard, test_data=test_om_hard, mtry=k)

      # k is only evaluated for this j if the number of soft clusters is reached
      if (j < nSoftClusts && soft_cluster_data$converged) {
        mse.cv.om_trate_soft_rf[i,j,k] <- train_mse_rf(train_data=train_om_soft, test_data=test_om_soft, mtry=k)
      }
    }
  }
}


# OM - SLOG 
# go through folds in repeat
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx

  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
  clusterward_hard <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")

  # hard coded for OM-Trate
  for (j in 2:nClusts) {
    
    # Create hard clusters
    hard_cluster_data <- hard_cluster_onehot(
      clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, 
      num_month_em_last_year = num_month_em_last_year, train_idx=train_idx, 
      test_idx=test_idx, dist_matrix = dists[[3]])

    train_om_hard <- hard_cluster_data$train_data
    test_om_hard <- hard_cluster_data$test_data
  
    # Create soft clusters
    soft_cluster_data <- soft_cluster(dist_matrix = dists[[3]], 
      train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
      num_month_em_last_year = num_month_em_last_year)

    train_om_soft <- soft_cluster_data$train_data
    test_om_soft <- soft_cluster_data$test_data

    for (k in 1:(floor(sqrt(nCovars + j - 2)))) {
     
      mse.cv.om_slog_hard_rf[i,j,k] <- train_mse_rf(train_data=train_om_hard, test_data=test_om_hard, mtry=k)
      
      # make sure we don't look at 25 columns (but we still want to look at mtry up to the number of soft clusters)
      # because it has multiple columns (one for each cluster)
      if (j < nSoftClusts && soft_cluster_data$converged) {
        mse.cv.om_slog_soft_rf[i,j,k] <- train_mse_rf(train_data=train_om_soft, test_data=test_om_soft, mtry=k)
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
    
    cluster_data <- hard_cluster_onehot(clusterward = clusterward_hard, nClusts = j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year, train_idx = train_idx, test_idx = test_idx, dist_matrix = dists[[2]])
    
    train_lcs_hard <- cluster_data$train_data
    test_lcs_hard <- cluster_data$test_data

    soft_cluster_data <- soft_cluster(dist_matrix = dists[[2]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, num_month_em_last_year = num_month_em_last_year)

    train_lcs_soft <- soft_cluster_data$train_data
    test_lcs_soft <- soft_cluster_data$test_data
    
    for (k in 1:(floor(sqrt(nCovars + j - 2)))) {

      mse.cv.lcs_hard_rf[i,j,k] <- train_mse_rf(train_data=train_lcs_hard, test_data=test_lcs_hard, mtry=k)

      if (j < nSoftClusts && soft_cluster_data$converged) {
        mse.cv.lcs_soft_rf[i,j,k] <- train_mse_rf(train_data=train_lcs_soft, test_data=test_lcs_soft, mtry=k)
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
#nSeqPcs <- ncol(rmetrics)

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
  
  for (j in 1:nSeqPcs) {
    train_rmets <- cbind(num_month_em_last_year[train_idx], mvad_covars[train_idx,], train_scores[,1:j,drop = FALSE])
    test_rmets <- cbind(num_month_em_last_year[test_idx], mvad_covars[test_idx,], test_scores[,1:j,drop = FALSE])
  
    colnames(train_rmets)[1] <- "num_month_em_last_year"
    colnames(test_rmets)[1] <- "num_month_em_last_year"
    
    for (k in 1:(floor(sqrt(nCovars + j - 1)))) {
      mse.cv.rmets_rf[i,j,k] <- train_mse_rf(train_data=train_rmets, test_data=test_rmets, mtry=k)
    }
  }
}

# Find the minimum performance for each fold-cluster combination 
# essentially this gets rid of all non-optimal performances (as meausured by mtry)
seq_mets_mins <- apply(mse.cv.rmets_rf, c(1, 2), min, na.rm = TRUE)
seq_mets_min_mtry <- apply(mse.cv.rmets_rf, c(1, 2), which.min)

# find average of best performances (by mtry) across folds 
seq_mets_pc_means_lm <- sqrt(apply(seq_mets_mins, 2, mean))
best_n_seq_mets_pc <- which.min(seq_mets_pc_means_lm)

# Old way for handling best way to do it 
# data frame to hold best combinations - long makes more sense 
mse.seq_clusts <- list()
mse.seq_clusts$om_trate_hard <- array(NA,c(folds,nClusts,nCovars + 1))
mse.seq_clusts$om_trate_soft <-  array(NA,c(folds,nClusts,nCovars + nSoftClusts - 1))
mse.seq_clusts$om_slog_hard <- array(NA,c(folds,nClusts,nCovars + 1))
mse.seq_clusts$om_slog_soft <- array(NA,c(folds,nClusts,nCovars + nSoftClusts - 1))
mse.seq_clusts$lcs_hard <- array(NA,c(folds,nClusts,nCovars + 1))
mse.seq_clusts$lcs_soft <- array(NA,c(folds,nClusts,nCovars + nSoftClusts - 1))


for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx
  
  # using rmetrics from above
  train_rmetrics <- rmetrics[train_idx, ]
  test_rmetrics <- rmetrics[test_idx, ]
  
  pca_comps_train <- prcomp(x=train_rmetrics, center=TRUE, scale=TRUE)
  
  train_scores <- pca_comps_train$x[,1:best_n_seq_mets_pc]
  test_scores <- predict(pca_comps_train, newdata = test_rmetrics)[,1:best_n_seq_mets_pc]
      
  # Do 6 different combinations of clustering and sequence PCS
  # create hierarchical clustering
  clusterward_hard_trate <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")
  clusterward_hard_lcs <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")
  clusterward_hard_slog <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")

  # create clusterings 


  for (j in 2:nClusts) {

    om_trate_clustered_hard <- hard_cluster_onehot(clusterward = clusterward_hard_trate, nClusts = j, 
      covars=mvad_covars, num_month_em_last_year, train_idx, test_idx, dists[[1]])
    om_trate_hard_train <- cbind(om_trate_clustered_hard$train_data, train_scores) 
    om_trate_hard_test <- cbind(om_trate_clustered_hard$test_data, test_scores) 
    
    # om-slog hard
  om_slog_clustered_hard <- hard_cluster_onehot(clusterward = clusterward_hard_slog, nClusts = j, 
      covars=mvad_covars, num_month_em_last_year, train_idx, test_idx, dists[[3]])
  om_slog_hard_train <- cbind(om_slog_clustered_hard$train_data, train_scores) 
  om_slog_hard_test <- cbind(om_slog_clustered_hard$test_data, test_scores)
    
  # lcs hard
  lcs_clustered_hard <- hard_cluster_onehot(clusterward = clusterward_hard_lcs, nClusts = j, 
      covars=mvad_covars, num_month_em_last_year, train_idx, test_idx, dists[[2]])
  lcs_hard_train <- cbind(lcs_clustered_hard$train_data, train_scores)
  lcs_hard_test <- cbind(lcs_clustered_hard$test_data, test_scores) 
    
  # om-trate soft
    if (j < nSoftClusts) {
      om_trate_clustered_soft <- soft_cluster(dists[[1]], train_idx, test_idx = test_idx,
        nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, num_month_em_last_year = num_month_em_last_year)
      om_trate_soft_train <- cbind(om_trate_clustered_soft$train_data, train_scores)
      om_trate_soft_test <- cbind(om_trate_clustered_soft$test_data, test_scores) 

      # om-slog soft
      om_slog_clustered_soft <- soft_cluster(dists[[3]], train_idx, test_idx = test_idx,
        nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, num_month_em_last_year = num_month_em_last_year)
      om_slog_soft_train <- cbind(om_slog_clustered_soft$train_data, train_scores)
      om_slog_soft_test <- cbind(om_slog_clustered_soft$test_data, test_scores)

      # lcs soft 
      lcs_clustered_soft <- soft_cluster(dists[[2]], train_idx, test_idx = test_idx,
        nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, num_month_em_last_year = num_month_em_last_year)
      lcs_soft_train <- cbind(lcs_clustered_soft$train_data, train_scores) 
      lcs_soft_test <- cbind(lcs_clustered_soft$test_data, test_scores) 
      
    }

    for (k in 1:(floor(sqrt(nCovars + j - 1)))) {

      mse.seq_clusts[["om_trate_hard"]][i, j,k] <- train_mse_rf(
        train_data = om_trate_hard_train, 
        test_data = om_trate_hard_test)
    
      mse.seq_clusts[["om_slog_hard"]][i, j,k] <- train_mse_rf(
        train_data = om_slog_hard_train, 
        test_data = om_slog_hard_test)

      mse.seq_clusts[["lcs_hard"]][i, j,k] <- train_mse_rf(
        train_data = lcs_hard_train, 
        test_data = lcs_hard_test)
      
      if (j < nSoftClusts) {
      
        if (om_trate_clustered_soft$converged) {
          mse.seq_clusts[["om_trate_soft"]][i,j,k] <- train_mse_rf(
            train_data = om_trate_soft_train, 
            test_data = om_trate_soft_test)
        }
      
        if (om_slog_clustered_soft$converged) {
          mse.seq_clusts[["om_slog_soft"]][i,j,k] <- train_mse_rf(
            train_data = om_slog_soft_train, 
            test_data = om_slog_soft_test)

        }

          if (lcs_clustered_soft$converged) {
          mse.seq_clusts[["lcs_soft"]][i,j,k] <- train_mse_rf(
            train_data = lcs_soft_train, 
            test_data = lcs_soft_test)
          }
        }
    }
  }
}


## Tracking for Sequence Metrics + Clustering 

# Get optimal performance for each seq pc + clustering method. 
min_mtry.seq_clusts <- list()
min_mtry.seq_clusts$om_trate_hard <- apply(mse.seq_clusts[[ "om_trate_hard"]], c(1,2), safe_which_min)
min_mtry.seq_clusts$om_trate_soft <-  apply(mse.seq_clusts[["om_trate_soft"]], c(1,2), safe_which_min)
min_mtry.seq_clusts$om_slog_hard <- apply(mse.seq_clusts[["om_slog_hard"]], c(1,2), safe_which_min)
min_mtry.seq_clusts$om_slog_soft <- apply(mse.seq_clusts[["om_slog_soft" ]], c(1,2), safe_which_min)
min_mtry.seq_clusts$lcs_hard <- apply(mse.seq_clusts[["lcs_hard"]], c(1,2), safe_which_min)
min_mtry.seq_clusts$lcs_soft <- apply(mse.seq_clusts[["lcs_soft"]], c(1,2), safe_which_min)

# Inf for n=1 clusters - ignore for now
min_mse.seq_clusts <- list()
min_mse.seq_clusts$om_trate_hard <- apply(mse.seq_clusts[[ "om_trate_hard"]], c(1,2), safe_min, na.rm = TRUE)
min_mse.seq_clusts$om_trate_soft <-  apply(mse.seq_clusts[["om_trate_soft"]], c(1,2),  safe_min, na.rm = TRUE)
min_mse.seq_clusts$om_slog_hard <- apply(mse.seq_clusts[["om_slog_hard"]], c(1,2),  safe_min, na.rm = TRUE)
min_mse.seq_clusts$om_slog_soft <- apply(mse.seq_clusts[["om_slog_soft" ]], c(1,2),  safe_min, na.rm = TRUE)
min_mse.seq_clusts$lcs_hard <- apply(mse.seq_clusts[["lcs_hard"]], c(1,2), safe_min, na.rm = TRUE)
min_mse.seq_clusts$lcs_soft <- apply(mse.seq_clusts[["lcs_soft"]], c(1,2), safe_min, na.rm = TRUE)

rmse.seqs_clusts <- list()
rmse.seqs_clusts$om_trate_hard <- sqrt(apply(min_mse.seq_clusts$om_trate_hard, 2, mean))
rmse.seqs_clusts$om_trate_soft <- sqrt(apply(min_mse.seq_clusts$om_trate_soft , 2, mean))
rmse.seqs_clusts$om_slog_hard <- sqrt(apply(min_mse.seq_clusts$om_slog_hard, 2, mean))
rmse.seqs_clusts$om_slog_soft <- sqrt(apply(min_mse.seq_clusts$om_slog_soft, 2, mean))
rmse.seqs_clusts$lcs_hard <- sqrt(apply(min_mse.seq_clusts$lcs_hard, 2, mean))
rmse.seqs_clusts$lcs_soft <- sqrt(apply(min_mse.seq_clusts$lcs_soft, 2, mean))

saveRDS(min_mtry.seq_clusts, "NonLinearSEQMBestMtry.rds")
saveRDS(min_mse.seq_clusts, "NonLinearSEQMMinMSE.rds")
saveRDS(rmse.seqs_clusts, "NonLinearSEQMAverRMSE.rds")


# Overall Performance Tracking 

# In our summary we want the average RMSE across the folds for each method for each number of clusters
# in other words we want to find the minimum performances (across the mtry) and then report the average of
# these




mse.comp <- list()
mse.comp$om_trate_hard <- mse.cv.om_trate_hard_rf
mse.comp$om_trate_soft <- mse.cv.om_trate_soft_rf
mse.comp$om_slog_hard  <- mse.cv.om_slog_hard_rf
mse.comp$om_slog_soft  <- mse.cv.om_slog_soft_rf
mse.comp$lcs_hard <- mse.cv.lcs_hard_rf
mse.comp$lcs_soft <- mse.cv.lcs_soft_rf
mse.comp$windows <- mse.cv.windows_rf
mse.comp$harm <- mse.cv.harm_rf



# keep actual mtry values to know optimal performances across fold and number of components 
min_mtry <- list()
min_mtry$harm_rf_min_mtry <- apply(mse.cv.harm_rf, c(1,2), which.min)
min_mtry$windows_rf_min_mtry <- apply(mse.cv.windows_rf , c(1,2), which.min)
min_mtry$om_trate_hard_rf_min_mtry <- apply(mse.cv.om_trate_hard_rf , c(1,2), which.min)
min_mtry$om_trate_soft_rf_min_mtry <- apply(mse.cv.om_trate_soft_rf,  c(1,2), which.min)
min_mtry$om_slog_hard_rf_min_mtry <- apply(mse.cv.om_slog_hard_rf , c(1,2), which.min)
min_mtry$om_slog_soft_rf_min_mtry <- apply(mse.cv.om_slog_soft_rf , c(1,2), which.min)
min_mtry$cs_hard_rf_min_mtry <- apply(mse.cv.lcs_hard_rf, c(1,2), which.min)
min_mtry$lcs_soft_rf_min_mtry <- apply(mse.cv.lcs_soft_rf, c(1,2), which.min)

# First find the minimum across dim 1 and 2 (this essentially gets the best performance
# by each combination)
min_mse <- list()
min_mse$harm <- apply(mse.cv.harm_rf, c(1,2), safe_min, na.rm = TRUE)
min_mse$windows <- apply(mse.cv.windows_rf , c(1,2), safe_min, na.rm = TRUE)
min_mse$om_trate_hard <- apply(mse.cv.om_trate_hard_rf , c(1,2), safe_min, na.rm = TRUE)
min_mse$om_trate_soft <- apply(mse.cv.om_trate_soft_rf,  c(1,2), safe_min, na.rm = TRUE)
min_mse$om_slog_hard <- apply(mse.cv.om_slog_hard_rf , c(1,2), safe_min, na.rm = TRUE)
min_mse$om_slog_soft <- apply(mse.cv.om_slog_soft_rf , c(1,2), safe_min, na.rm = TRUE)
min_mse$lcs_hard <- apply(mse.cv.lcs_hard_rf, c(1,2), safe_min, na.rm = TRUE)
min_mse$lcs_soft <- apply(mse.cv.lcs_soft_rf, c(1,2), safe_min, na.rm = TRUE)

# get means across the folds (of the best performances by mtry)
rmse.comp <- list()
rmse.comp$harm <- sqrt(apply(min_mse$harm, 2, mean))
rmse.comp$windows <- sqrt(apply(min_mse$windows, 2, mean))
rmse.comp$om_trate_hard <- sqrt(apply(min_mse$om_trate_hard, 2, mean))
rmse.comp$om_trate_soft <- sqrt(apply(min_mse$om_trate_soft, 2, mean))
rmse.comp$om_slog_hard <- sqrt(apply(min_mse$om_slog_hard, 2, mean))
rmse.comp$om_slog_soft <- sqrt(apply(min_mse$om_slog_soft, 2, mean))
rmse.comp$lcs_hard <- sqrt(apply(min_mse$lcs_hard, 2, mean))
rmse.comp$lcs_soft <- sqrt(apply(min_mse$lcs_soft, 2, mean))




saveRDS(mse.comp, "NonLinearCompFullMSEs.rds")
saveRDS(min_mtry, "NonLinearCompBestMtry.rds")
saveRDS(min_mse, "NonLinearCompBestMSE.rds")
saveRDS(rmse.comp, "NonLinearCompAverRMSES.rds")


method_labels <- c("om_trate_hard"="OM T-Rate (Hard)", 
"om_trate_soft"="OM T-Rate (Soft)", 
"om_slog_hard"="OM INDELSLOG (Hard)", 
"om_slog_soft"="OM INDELSLOG (Soft)",
"lcs_hard"="LCS (Hard)", 
"lcs_soft"="LCS (Soft)", 
"windows"="Windows", 
"harm"="CFDA")

rmse_plot_data <- imap_dfr(rmse.comp, ~ data.frame(
  index = seq_along(.x),
  value = .x,
  group = .y
))
colnames(rmse_plot_data) <- c("comps", "rmse", "method")
rmse_plot_data$method <- method_labels[rmse_plot_data$method]


pdf("NonLinearCompPlot.pdf",width=8,height=6) 
ggplot(data = rmse_plot_data, 
       aes(x = comps, y = rmse, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  labs(title = "Random Forest Performance by Method",
       x = "Components/Clusters",
       y = "RMSE Value (Averaged)",
       color = "Method") +
  theme_minimal()

dev.off()

saveRDS(rmse_plot_data, file="NonLinearCompetitionRMSE.rds")



  







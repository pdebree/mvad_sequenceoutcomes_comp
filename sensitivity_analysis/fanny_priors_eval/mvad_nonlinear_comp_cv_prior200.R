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
library(foreach)
library(doParallel)

source("seqout_utils.R")


# Detect cores allocated by Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# Register the cluster
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 123) # for reproducability
registerDoParallel(cl)


folds <- 5
nClusts <- 25 
nWindows <- 16 
nHarms <- 25
nSeqPcs <- 12
fuzz_soft <- 1.5
nSoftClusts <- 13 
nCovars <- 11 
prior.n.equiv <- 200

nCVs <- 20

comp <- list()
comp$om_trate_hard <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$om_trate_soft <- array(NA,c(folds,nSoftClusts,nCovars + nSoftClusts - 2))
comp$om_slog_hard <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$om_slog_soft <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
comp$lcs_hard <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$lcs_soft <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
comp$windows <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))
comp$harm <- array(NA,c(folds,nHarms, nCovars + nHarms - 1 ))
comp$demos <- array(NA, c(folds, nCovars - 1))

cv_comp <- list()
cv_comp$mse <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$cov <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$mpiw <- replicate(nCVs, comp, simplify = FALSE)

task_vec <- data.frame(cv_index = 1:nCVs)

results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras", "fpc", "NbClust", "ranger", "tuneRanger", "gt")) %dopar% {
  
  cv_index <- task_vec$cv_index[m]

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
  prior.n.equiv <- 200


  # Arrays for holding outcomes fits 
  # soft clusters - 2 (because we never look at index=1 clusters and indexing is always from 1)
  mpiw.cv.harm_rf <- cov.cv.harm_rf <- mse.cv.harm_rf <- array(NA,c(folds,nHarms, nCovars + nHarms - 1 ))
  mpiw.cv.windows_rf <- cov.cv.windows_rf <- mse.cv.windows_rf <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))
  mpiw.cv.om_trate_hard_rf <- cov.cv.om_trate_hard_rf <- mse.cv.om_trate_hard_rf <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
  mpiw.cv.om_trate_soft_rf <- cov.cv.om_trate_soft_rf <- mse.cv.om_trate_soft_rf <- array(NA,c(folds,nSoftClusts,nCovars + nSoftClusts - 2))
  mpiw.cv.om_slog_hard_rf <- cov.cv.om_slog_hard_rf <- mse.cv.om_slog_hard_rf <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
  mpiw.cv.om_slog_soft_rf <- cov.cv.om_slog_soft_rf <- mse.cv.om_slog_soft_rf <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
  mpiw.cv.lcs_hard_rf <- cov.cv.lcs_hard_rf <- mse.cv.lcs_hard_rf <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
  mpiw.cv.lcs_soft_rf <- cov.cv.lcs_soft_rf <- mse.cv.lcs_soft_rf <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
  mpiw.cv.rmets_rf <- cov.cv.rmets_rf  <- mse.cv.rmets_rf <- array(NA,c(folds, nSeqPcs, nCovars + nSeqPcs - 1))
  mpiw.cv.demos <- cov.cv.demos  <- mse.cv.demos <- array(NA,c(folds, nCovars - 1))



  ## OM-TRATE
  # go through folds in repeat
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx

    # hard coded for OM-Trate
    for (j in 2:nClusts) {
    
      # Create soft clusters
      soft_cluster_data <- soft_cluster(dist_matrix = dists[[1]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, y = num_month_em_last_year, prior.n.equiv=prior.n.equiv)

      train_om_soft <- soft_cluster_data$train_data
      test_om_soft <- soft_cluster_data$test_data

      # have to minus 2 because of indexing j from 2 (instead of 1 - can't look at 1 soft cluster)
      for (k in 1:(nCovars + j - 2)) {
        # k is only evaluated for this j if the number of soft clusters is reached
        if (j < nSoftClusts && soft_cluster_data$converged) {
          soft_om_fit <- fit_rf(train_data=train_om_soft, test_data=test_om_soft, mtry=k)
          mse.cv.om_trate_soft_rf[i,j,k] <- soft_om_fit$mse
          cov.cv.om_trate_soft_rf[i, j, k] <- soft_om_fit$cov
          mpiw.cv.om_trate_soft_rf[i, j, k] <- soft_om_fit$mpiw

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
      
    
      # Create soft clusters
      soft_cluster_data <- soft_cluster(dist_matrix = dists[[3]], 
        train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, prior.n.equiv=prior.n.equiv)

      train_om_soft <- soft_cluster_data$train_data
      test_om_soft <- soft_cluster_data$test_data

      for (k in 1:(nCovars + j - 2)) {

        if (j < nSoftClusts && soft_cluster_data$converged) {
          soft_slog_fit <- fit_rf(train_data=train_om_soft, test_data=test_om_soft, mtry=k)
          mse.cv.om_slog_soft_rf[i,j,k] <- soft_slog_fit$mse
          cov.cv.om_slog_soft_rf[i, j, k] <- soft_slog_fit$cov
          mpiw.cv.om_slog_soft_rf[i, j, k] <- soft_slog_fit$mpiw
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


    for (j in 2:nClusts) {
      

      soft_cluster_data <- soft_cluster(dist_matrix = dists[[2]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, y = num_month_em_last_year, prior.n.equiv=prior.n.equiv)

      train_lcs_soft <- soft_cluster_data$train_data
      test_lcs_soft <- soft_cluster_data$test_data
      
      for (k in 1:(nCovars + j - 2)) {


        if (j < nSoftClusts && soft_cluster_data$converged) {
          soft_lcs_fit <- fit_rf(train_data=train_lcs_soft, test_data=test_lcs_soft, mtry=k)
          mse.cv.lcs_soft_rf[i,j,k] <- soft_lcs_fit$mse
          cov.cv.lcs_soft_rf[i, j, k] <- soft_lcs_fit$cov
          mpiw.cv.lcs_soft_rf[i, j, k] <- soft_lcs_fit$mpiw
        }
      }
    }
  }
    

  mse.comp <- list()
  mse.comp$om_trate_hard <-mse.cv.om_trate_hard_rf  
  mse.comp$om_trate_soft <- mse.cv.om_trate_soft_rf 
  mse.comp$om_slog_hard  <- mse.cv.om_slog_hard_rf 
  mse.comp$om_slog_soft  <- mse.cv.om_slog_soft_rf 
  mse.comp$lcs_hard <- mse.cv.lcs_hard_rf 
  mse.comp$lcs_soft <- mse.cv.lcs_soft_rf 
  mse.comp$windows <- mse.cv.windows_rf 
  mse.comp$harm <- mse.cv.harm_rf 
  mse.comp$demos <- mse.cv.demos

  
  cov.comp <- list()
  cov.comp$om_trate_hard <- cov.cv.om_trate_hard_rf 
  cov.comp$om_trate_soft <- cov.cv.om_trate_soft_rf 
  cov.comp$om_slog_hard  <- cov.cv.om_slog_hard_rf 
  cov.comp$om_slog_soft  <- cov.cv.om_slog_soft_rf 
  cov.comp$lcs_hard <- cov.cv.lcs_hard_rf 
  cov.comp$lcs_soft <- cov.cv.lcs_soft_rf 
  cov.comp$windows <- cov.cv.windows_rf 
  cov.comp$harm <- cov.cv.harm_rf
  cov.comp$demos <- cov.cv.demos

  mpiw.comp <- list()
  mpiw.comp$om_trate_hard <-  mpiw.cv.om_trate_hard_rf 
  mpiw.comp$om_trate_soft <-  mpiw.cv.om_trate_soft_rf 
  mpiw.comp$om_slog_hard  <- mpiw.cv.om_slog_hard_rf 
  mpiw.comp$om_slog_soft  <-  mpiw.cv.om_slog_soft_rf 
  mpiw.comp$lcs_hard <- mpiw.cv.lcs_hard_rf 
  mpiw.comp$lcs_soft <- mpiw.cv.lcs_soft_rf 
  mpiw.comp$windows <- mpiw.cv.windows_rf 
  mpiw.comp$harm <-  mpiw.cv.harm_rf 
  mpiw.comp$demos <- mpiw.cv.demos

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
    
}

for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_NonLinearComp_prior200.rds")

  
# Stop the cluster
stopCluster(cl)

  
  
  
  
  
  
  
  
  
  
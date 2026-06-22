# CFDA Linear Competition 
# First written for Multi-Level Models practicum in Fall 2025

library(tidyverse)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)

source("seqout_utils.R")

# # Detect cores allocated by Slurm
# n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# # Register the cluster
# cl <- makeCluster(n_cores)
# clusterSetRNGStream(cl, iseed = 123) # for reproducability
# registerDoParallel(cl)


folds <- 5
nClusts <- 25 
nWindows <- 16 
nHarms <- 25
nSeqPcs <- 12
fuzz_soft <- 1.5
nSoftClusts <- 13 
prior.n.equiv <- 200

nCVs <- 20

comp <- list()
comp$om_trate_hard <- array(NA,c(folds,nClusts))
comp$om_trate_soft <- array(NA,c(folds,nClusts))
comp$om_slog_hard <- array(NA,c(folds,nClusts))
comp$om_slog_soft <- array(NA,c(folds,nClusts))
comp$lcs_hard <- array(NA,c(folds,nClusts))
comp$lcs_soft <- array(NA,c(folds,nClusts))
comp$windows <- array(NA,c(folds,nWindows))
comp$harm <- array(NA,c(folds,nHarms))
comp$demos <- array(NA, c(folds))

cv_comp <- list()
cv_comp$mse <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$cov <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$mpiw <- replicate(nCVs, comp, simplify = FALSE)


task_vec <- data.frame(cv_index = 1:nCVs)


results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {

  cv_index <- task_vec$cv_index[m]
  
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

  ### Training set up 
  folds <-  5 
  nrecs <- nrow(mvad[1])

  # Set up training and test split so it is the same across the competition 
  idx.orig <- rep(1:folds,each=floor(nrecs/folds))
  if (nrecs %% folds != 0) idx.orig <- c(idx.orig,1:(nrecs %% folds))
  idx <- sample(idx.orig) #shuffle

  # Creates a random shuffle of the indices to be used in the the train-test split and folds. 
  shuffle_index <- sample(1:nrecs, nrecs)
  ids <- sort(unique(mvad$id))


  # Number of each method to check
  nClusts <- 25 
  nWindows <- 16 
  nHarms <- 25
  nSeqPcs <- 12

  fuzz_soft <- 1.5
  nSoftClusts <- 13 

  ### Arrays to hold RMSEs from different runs 
  mse.cv.om_trate_hard <- cov.cv.om_trate_hard <- mpiw.cv.om_trate_hard <- array(NA,c(folds,nClusts))
  mse.cv.om_trate_soft <- cov.cv.om_trate_soft <- mpiw.cv.om_trate_soft <- array(NA,c(folds,nClusts))
  mse.cv.om_slog_hard <- cov.cv.om_slog_hard <-  mpiw.cv.om_slog_hard <- array(NA,c(folds,nClusts))
  mse.cv.om_slog_soft <- cov.cv.om_slog_soft <- mpiw.cv.om_slog_soft <- array(NA,c(folds,nClusts))
  mse.cv.lcs_hard <- cov.cv.lcs_hard <- mpiw.cv.lcs_hard <- array(NA,c(folds,nClusts))
  mse.cv.lcs_soft <- cov.cv.lcs_soft <- mpiw.cv.lcs_soft <- array(NA,c(folds,nClusts))
  mse.cv.windows <- cov.cv.windows <- mpiw.cv.windows <- array(NA,c(folds,nWindows))
  mse.cv.harm <- cov.cv.harm  <- mpiw.cv.harm <- array(NA,c(folds,nHarms))
  mse.cv.mets <- cov.cv.mets <- mpiw.cv.mets  <- array(NA,c(folds,nSeqPcs))
  mse.cv.demos <- cov.cv.demos <- mpiw.cv.demos  <- array(NA,c(folds))



  ### Cross Validation for Clustering Methods with OM-Transition Rate 
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx
    

    
    # hard coded for OM-Trate
    for (j in 2:nClusts) {
      
      #soft clustering - fanny 
      if (j < nSoftClusts) {
        soft_cluster_data <- soft_cluster(dist_matrix = dists[[1]], 
          train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
          y = num_month_em_last_year, prior.n.equiv=prior.n.equiv)

        soft_output <- fit_linear(soft_cluster_data$train_data, soft_cluster_data$test_data)
        
        if (soft_cluster_data$converged) {
          mse.cv.om_trate_soft[i,j] <- soft_output$mse
          cov.cv.om_trate_soft[i,j] <- soft_output$coverage
          mpiw.cv.om_trate_soft[i,j] <- soft_output$mpiw

          
      
        }
      }
    }
  }

  ### Cross Validation for Clustering Methods with OM-Slog Rate 
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx
    

    
    # hard coded for OM-Trate
    for (j in 2:nClusts) {
      
      #soft clustering - fanny 
      if (j < nSoftClusts) {
        soft_cluster_data <- soft_cluster(dist_matrix = dists[[3]], 
          train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
          y = num_month_em_last_year, prior.n.equiv=prior.n.equiv)

        soft_output <- fit_linear(soft_cluster_data$train_data, soft_cluster_data$test_data)
        
        if (soft_cluster_data$converged) {
          mse.cv.om_slog_soft[i,j] <- soft_output$mse
          cov.cv.om_slog_soft[i,j] <- soft_output$coverage
          mpiw.cv.om_slog_soft[i,j] <- soft_output$mpiw
        }
      }
    }
  }


  ### Cross Validation for Clustering Methods with LCS as a Distance Measure
  for (i in 1:folds) {
    test_idx <- idx == i
    train_idx <- !test_idx
    
    
    for (j in 2:nClusts) {
      
      #soft clustering - fanny 
      if (j < nSoftClusts) {
        soft_cluster_data <- soft_cluster(dist_matrix = dists[[2]], 
          train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
          y = num_month_em_last_year, prior.n.equiv=prior.n.equiv)

        soft_output <- fit_linear(soft_cluster_data$train_data, soft_cluster_data$test_data)
        
        if (soft_cluster_data$converged) {
          mse.cv.lcs_soft[i,j] <- soft_output$mse
          cov.cv.lcs_soft[i,j] <- soft_output$coverage
          mpiw.cv.lcs_soft[i,j] <- soft_output$mpiw
        }
      }
    }
  }


  mse.comp <- list()
  mse.comp$om_trate_hard <- mse.cv.om_trate_hard 
  mse.comp$om_trate_soft <- mse.cv.om_trate_soft 
  mse.comp$om_slog_hard  <- mse.cv.om_slog_hard 
  mse.comp$om_slog_soft  <- mse.cv.om_slog_soft 
  mse.comp$lcs_hard <- mse.cv.lcs_hard
  mse.comp$lcs_soft <- mse.cv.lcs_soft
  mse.comp$windows <- mse.cv.windows
  mse.comp$harm <- mse.cv.harm
  mse.comp$mets <- mse.cv.mets 
  mse.comp$demos <- mse.cv.demos

  cov.comp <- list()
  cov.comp$om_trate_hard <- cov.cv.om_trate_hard 
  cov.comp$om_trate_soft <- cov.cv.om_trate_soft 
  cov.comp$om_slog_hard  <- cov.cv.om_slog_hard 
  cov.comp$om_slog_soft  <- cov.cv.om_slog_soft 
  cov.comp$lcs_hard <- cov.cv.lcs_hard
  cov.comp$lcs_soft <- cov.cv.lcs_soft
  cov.comp$windows <- cov.cv.windows
  cov.comp$harm <- cov.cv.harm
  cov.comp$mets <- cov.cv.mets 
  cov.comp$demos <- cov.cv.demos


  mpiw.comp <- list()
  mpiw.comp$om_trate_hard <- mpiw.cv.om_trate_hard 
  mpiw.comp$om_trate_soft <- mpiw.cv.om_trate_soft 
  mpiw.comp$om_slog_hard  <- mpiw.cv.om_slog_hard 
  mpiw.comp$om_slog_soft  <- mpiw.cv.om_slog_soft 
  mpiw.comp$lcs_hard <- mpiw.cv.lcs_hard
  mpiw.comp$lcs_soft <- mpiw.cv.lcs_soft
  mpiw.comp$windows <- mpiw.cv.windows
  mpiw.comp$harm <- mpiw.cv.harm
  mpiw.comp$mets <- mpiw.cv.mets 
  mpiw.comp$demos <- mpiw.cv.demos

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
}

# Populate global matrices from return list 



for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_LinearComp_200prior.rds")

  
# Stop the cluster
stopCluster(cl)



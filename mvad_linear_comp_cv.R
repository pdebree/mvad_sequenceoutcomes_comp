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


# Detect cores for remote running 
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 123) # for reproducability
registerDoParallel(cl)

folds <- 5
nCVs <- 20

nClusts <- 25 
nWindows <- 16 
nHarms <- 25
nCovars <- 11 
nSeqMetsPCs <- 3

fuzz_soft <- 1.75
max_soft_iters <- 2000
tolerance <- 1e-10
knn_soft_assign <- 15


# Create data forms to hold CV outputs 
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
comp$windows_3pc <- array(NA,c(folds,nWindows))
comp$harm_3pc <- array(NA,c(folds,nHarms))
comp$om_trate_soft_3pc <- array(NA,c(folds, nClusts))
comp$om_slog_soft_3pc <- array(NA,c(folds, nClusts))
comp$lcs_soft_3pc <- array(NA,c(folds, nClusts))
comp$om_slog_hard_3pc <- array(NA,c(folds,nClusts))

con.comp <- list()
con.comp$lcs_soft <- array(NA, c(folds, nClusts))
con.comp$om_trate_soft <- array(NA, c(folds, nClusts))
con.comp$om_slog_soft <- array(NA, c(folds, nClusts))

cv_comp <- list()
cv_comp$mse <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$cov <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$mpiw <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$convergence <- replicate(nCVs, con.comp, simplify = FALSE)

task_vec <- data.frame(cv_index = 1:nCVs)


# Iterate over tasks for each 
results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {

  cv_index <- task_vec$cv_index[m]
  
  ### Data Preparation 
  data(mvad)
  mvad_last_year <- mvad[,75:86]
  num_month_em_last_year <- apply(mvad_last_year, 1, function(x) length(which(x=="employment")))
  mvad_covars <- mvad[3:14] %>% dplyr::select(-Western) #reference group

  # Create distance matrices 
  mvad.seq <- create_seq_data(mvad)
  dists <- create_dists(mvad.seq)

  # Set up training and test split so it is the same across the cross validation
  nrecs <- nrow(mvad[1])
  idx.orig <- rep(1:folds,each=floor(nrecs/folds))
  if (nrecs %% folds != 0) idx.orig <- c(idx.orig,1:(nrecs %% folds))
  idx <- sample(idx.orig) #shuffle


  ### Arrays to hold RMSEs from different runs 
  mse.cv.om_trate_hard <- cov.cv.om_trate_hard <- mpiw.cv.om_trate_hard <- array(NA,c(folds,nClusts))
  mse.cv.om_trate_soft <- cov.cv.om_trate_soft <- mpiw.cv.om_trate_soft <- array(NA,c(folds,nClusts))
  mse.cv.om_slog_hard <- cov.cv.om_slog_hard <-  mpiw.cv.om_slog_hard <- array(NA,c(folds,nClusts))
  mse.cv.om_slog_soft <- cov.cv.om_slog_soft <- mpiw.cv.om_slog_soft <- array(NA,c(folds,nClusts))
  mse.cv.lcs_hard <- cov.cv.lcs_hard <- mpiw.cv.lcs_hard <- array(NA,c(folds,nClusts))
  mse.cv.lcs_soft <- cov.cv.lcs_soft <- mpiw.cv.lcs_soft <- array(NA,c(folds,nClusts))
  mse.cv.windows <- cov.cv.windows <- mpiw.cv.windows <- array(NA,c(folds,nWindows))
  mse.cv.harm <- cov.cv.harm  <- mpiw.cv.harm <- array(NA,c(folds,nHarms))
  mse.cv.demos <- cov.cv.demos <- mpiw.cv.demos  <- array(NA,c(folds))

  # For adding in Sequence Metric Principal Components
  mpiw.cv.om_slog_hard_3pc <- cov.cv.om_slog_hard_3pc <- mse.cv.om_slog_hard_3pc  <- array(NA,c(folds,nClusts)) 
  mpiw.cv.om_trate_soft_3pc  <- cov.cv.om_trate_soft_3pc  <- mse.cv.om_trate_soft_3pc <- array(NA,c(folds, nClusts))
  mpiw.cv.om_slog_soft_3pc  <- cov.cv.om_slog_soft_3pc  <- mse.cv.om_slog_soft_3pc <- array(NA,c(folds, nClusts))
  mpiw.cv.lcs_soft_3pc  <- cov.cv.lcs_soft_3pc  <- mse.cv.lcs_soft_3pc <- array(NA,c(folds, nClusts))
  mse.cv.windows_3pc <- cov.cv.windows_3pc <- mpiw.cv.windows_3pc <- array(NA,c(folds,nWindows))
  mse.cv.harm_3pc <- cov.cv.harm_3pc  <- mpiw.cv.harm_3pc <- array(NA,c(folds,nHarms))

  # for tracking convergence / condition number 
  con.cv.lcs_soft <- con.cv.om_trate_soft <- con.cv.om_slog_soft <- array(NA, c(folds, nClusts))

  # Read data for different methods
  rmetrics <- create_rmetrics(mvad.seq=mvad.seq)
  year_state_counts <- create_year_state_counts(mvad)
  mvad_long <- create_long_data(mvad)
  basis <- create.bspline.basis(c(0, 35), nbasis = 6, norder = 4)

  ### Cross Validation for Demographics only
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx
    
    # Demos 
    demos_data <- mvad_covars |> mutate(y = num_month_em_last_year)
    demos_output <- fit_linear(demos_data[train_idx,], demos_data[test_idx,])
    mse.cv.demos[i] <- demos_output$mse
    cov.cv.demos[i] <- demos_output$coverage
    mpiw.cv.demos[i] <- demos_output$mpiw


    # Make principal components from sequence metrics pcs 
    seqmets_pcs <- create_seqmets_pcs(rmetrics, train_idx, test_idx, nSeqMetsPCs)
    train_seq_scores <- seqmets_pcs$train_pcs
    test_seq_scores <- seqmets_pcs$test_pcs
    
    # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM-Slog
    clusterward_hard <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")

    # hard coded for OM-SLOG
    for (j in 2:nClusts) {

      hard_cluster_data <- hard_cluster(
        clusterward=clusterward_hard, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, train_idx=train_idx, 
        test_idx=test_idx, dist_matrix = dists[[3]])

    
      hard_output <- fit_linear(hard_cluster_data$train_data, hard_cluster_data$test_data)

      # hard clustering
      mse.cv.om_slog_hard[i,j] <- hard_output$mse
      cov.cv.om_slog_hard[i,j] <- hard_output$coverage
      mpiw.cv.om_slog_hard[i,j] <- hard_output$mpiw


      # Sequnece Metrics calcualtions
      train_hard_seq <- cbind(hard_cluster_data$train_data, train_seq_scores)
      test_hard_seq <- cbind(hard_cluster_data$test_data, test_seq_scores)

      hard_seq_output <- fit_linear(train_hard_seq, test_hard_seq)

      mse.cv.om_slog_hard_3pc[i,j] <- hard_seq_output$mse
      cov.cv.om_slog_hard_3pc[i,j] <- hard_seq_output$coverage
      mpiw.cv.om_slog_hard_3pc[i,j] <- hard_seq_output$mpiw
        
      
      # soft clustering
      soft_cluster_data_om_slog <- soft_cluster(dist_matrix = dists[[3]], 
        train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, max_iters = max_soft_iters, fuzziness = fuzz_soft,
        knn_soft_assign = knn_soft_assign, tol=tolerance)
      
      con.cv.om_slog_soft[i, j] <- soft_cluster_data_om_slog$converged

      train_om_slog_soft_seq <- cbind(soft_cluster_data_om_slog$train_data, train_seq_scores)
      test_om_slog_soft_seq <- cbind(soft_cluster_data_om_slog$test_data, test_seq_scores)

      if (soft_cluster_data_om_slog$converged) {

        soft_om_slog_output <- fit_linear(soft_cluster_data_om_slog$train_data, soft_cluster_data_om_slog$test_data)

        mse.cv.om_slog_soft[i,j] <- soft_om_slog_output$mse
        cov.cv.om_slog_soft[i,j] <- soft_om_slog_output$coverage
        mpiw.cv.om_slog_soft[i,j] <- soft_om_slog_output$mpiw

        om_slog_soft_seq_output <- fit_linear(train_om_slog_soft_seq, test_om_slog_soft_seq)

        mse.cv.om_slog_soft_3pc[i,j] <- om_slog_soft_seq_output$mse
        cov.cv.om_slog_soft_3pc[i,j] <- om_slog_soft_seq_output$coverage
        mpiw.cv.om_slog_soft_3pc[i,j]  <- om_slog_soft_seq_output$mpiw 
      
      }
    }
    
    # Make PCs for fold 
    mvad_count_pcs <- create_counts_pcs(year_state_counts, mvad_covars, train_idx, test_idx, num_month_em_last_year)
    train_mvad_windows <- mvad_count_pcs$train_data
    test_mvad_windows <- mvad_count_pcs$test_data
  

    for (j in 1:nWindows) {

      train_j_wind <- train_mvad_windows[, 1:c(12+j)]
      test_j_wind <- test_mvad_windows[, 1:c(12+j)]

      wind_fit <- fit_linear(train_j_wind, test_j_wind)
      mse.cv.windows[i,j] <- wind_fit$mse
      cov.cv.windows[i,j] <- wind_fit$coverage
      mpiw.cv.windows[i,j] <- wind_fit$mpiw


      # Add in sequence metrics 
      wind_fit_seq <- fit_linear(cbind(train_j_wind, train_seq_scores), cbind(test_j_wind, test_seq_scores))
      mse.cv.windows_3pc[i,j] <- wind_fit_seq$mse
      cov.cv.windows_3pc[i,j] <- wind_fit_seq$coverage
      mpiw.cv.windows_3pc[i,j] <- wind_fit_seq$mpiw

    }
  
    # Create harmonics data for folds 
    harm_data <- create_cfda_harms(mvad_long, train_idx, test_idx, basis)
    train_harm <- harm_data$train_data
    test_harm <- harm_data$test_data


    for (j in 1:nHarms) {

      train_j_harm <- train_harm[, 1:c(12+j)]
      test_j_harm <- test_harm[, 1:c(12+j)]
      
      harm_fit <- fit_linear(train_j_harm, test_j_harm)
      
      mse.cv.harm[i,j] <- harm_fit$mse
      cov.cv.harm[i,j] <- harm_fit$coverage
      mpiw.cv.harm[i,j] <- harm_fit$mpiw
      
      # Add in Sequence metrics 
      harm_fit_seq <- fit_linear(cbind(train_j_harm, train_seq_scores), cbind(test_j_harm, test_seq_scores))
      mse.cv.harm_3pc[i,j] <- harm_fit_seq$mse
      cov.cv.harm_3pc[i,j] <- harm_fit_seq$coverage
      mpiw.cv.harm_3pc[i,j] <- harm_fit_seq$mpiw

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
  mse.comp$demos <- mse.cv.demos
  mse.comp$om_trate_soft_3pc <- mse.cv.om_trate_soft_3pc
  mse.comp$om_slog_soft_3pc <- mse.cv.om_slog_soft_3pc
  mse.comp$lcs_soft_3pc <- mse.cv.lcs_soft_3pc
  mse.comp$om_slog_hard_3pc <- mse.cv.om_slog_hard_3pc 
  mse.comp$harm_3pc <- mse.cv.harm_3pc
  mse.comp$windows_3pc <- mse.cv.windows_3pc
  
 

  cov.comp <- list()
  cov.comp$om_trate_hard <- cov.cv.om_trate_hard 
  cov.comp$om_trate_soft <- cov.cv.om_trate_soft 
  cov.comp$om_slog_hard  <- cov.cv.om_slog_hard 
  cov.comp$om_slog_soft  <- cov.cv.om_slog_soft 
  cov.comp$lcs_hard <- cov.cv.lcs_hard
  cov.comp$lcs_soft <- cov.cv.lcs_soft
  cov.comp$windows <- cov.cv.windows
  cov.comp$harm <- cov.cv.harm
  cov.comp$demos <- cov.cv.demos
  cov.comp$om_trate_soft_3pc <- cov.cv.om_trate_soft_3pc
  cov.comp$om_slog_soft_3pc <- cov.cv.om_slog_soft_3pc
  cov.comp$lcs_soft_3pc <- cov.cv.lcs_soft_3pc
  cov.comp$om_slog_hard_3pc <- cov.cv.om_slog_hard_3pc
  cov.comp$harm_3pc <- cov.cv.harm_3pc
  cov.comp$windows_3pc <- cov.cv.windows_3pc
  
  mpiw.comp <- list()
  mpiw.comp$om_trate_hard <- mpiw.cv.om_trate_hard 
  mpiw.comp$om_trate_soft <- mpiw.cv.om_trate_soft 
  mpiw.comp$om_slog_hard  <- mpiw.cv.om_slog_hard 
  mpiw.comp$om_slog_soft  <- mpiw.cv.om_slog_soft 
  mpiw.comp$lcs_hard <- mpiw.cv.lcs_hard
  mpiw.comp$lcs_soft <- mpiw.cv.lcs_soft
  mpiw.comp$windows <- mpiw.cv.windows
  mpiw.comp$harm <- mpiw.cv.harm
  mpiw.comp$demos <- mpiw.cv.demos
  mpiw.comp$om_trate_soft_3pc <- mpiw.cv.om_trate_soft_3pc  
  mpiw.comp$om_slog_soft_3pc <- mpiw.cv.om_slog_soft_3pc  
  mpiw.comp$lcs_soft_3pc <- mpiw.cv.lcs_soft_3pc  
  mpiw.comp$om_slog_hard_3pc <- mpiw.cv.om_slog_hard_3pc
  mpiw.comp$harm_3pc <- mpiw.cv.harm_3pc
  mpiw.comp$windows_3pc <- mpiw.cv.windows_3pc

  con.comp <- list()
  con.comp$lcs_soft <- con.cv.lcs_soft
  con.comp$om_trate_soft <- con.cv.om_trate_soft
  con.comp$om_slog_soft <- con.cv.om_slog_soft

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, con.comp=con.comp, cv_index=cv_index)
  }


for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
  cv_comp$convergence[[n]] <- result$con.comp
}


saveRDS(cv_comp, file="MVADLinearComp_Results.rds")

  
# Stop the cluster
stopCluster(cl)



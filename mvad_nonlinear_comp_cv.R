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

comp <- list()

# here mtry for hard clustering becomes nCovers + 1 - 1 
comp$om_trate_hard <- array(NA,c(folds,nClusts, nCovars))
comp$om_trate_soft <- array(NA,c(folds,nClusts,nCovars + nClusts - 2))
comp$om_slog_hard <- array(NA,c(folds,nClusts, nCovars))
comp$om_slog_soft <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$lcs_hard <- array(NA,c(folds,nClusts, nCovars))
comp$lcs_soft <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$windows <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))
comp$harm <- array(NA,c(folds,nHarms, nCovars + nHarms - 1 ))
comp$demos <- array(NA, c(folds, nCovars - 1))
comp$om_trate_soft_3pc <- array(NA,c(folds, nClusts, nCovars + nClusts + 3 - 2))
comp$om_slog_soft_3pc <- array(NA,c(folds, nClusts, nCovars + nClusts + 3 - 2))
comp$lcs_soft_3pc <- array(NA,c(folds, nClusts, nCovars + nClusts + 3 - 2))
comp$om_slog_hard_3pc <- array(NA,c(folds,nClusts, nCovars + 3))
comp$windows_3pc <- array(NA,c(folds,nWindows, nCovars + nWindows + 3 - 1))
comp$harm_3pc <- array(NA,c(folds,nHarms, nCovars + nHarms  + 3 - 1 ))


# Don't need to add in Mtry because we will do multiply mtry vals for the same fit
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

results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras", "fpc", "NbClust", "ranger", "tuneRanger", "gt")) %dopar% {
  
  cv_index <- task_vec$cv_index[m]

  # Data Preparation 
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


  # Arrays for holding outcomes fits 
  # soft clusters - 2 (because we never look at index=1 clusters and indexing is always from 1)
  mpiw.cv.harm_rf <- cov.cv.harm_rf <- mse.cv.harm_rf <- array(NA,c(folds,nHarms, nCovars + nHarms - 1 ))
  mpiw.cv.windows_rf <- cov.cv.windows_rf <- mse.cv.windows_rf <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))
  mpiw.cv.om_trate_hard_rf <- cov.cv.om_trate_hard_rf <- mse.cv.om_trate_hard_rf <- array(NA,c(folds,nClusts, nCovars))
  mpiw.cv.om_trate_soft_rf <- cov.cv.om_trate_soft_rf <- mse.cv.om_trate_soft_rf <- array(NA,c(folds,nClusts,nCovars + nClusts - 2))
  mpiw.cv.om_slog_hard_rf <- cov.cv.om_slog_hard_rf <- mse.cv.om_slog_hard_rf <- array(NA,c(folds,nClusts, nCovars))
  mpiw.cv.om_slog_soft_rf <- cov.cv.om_slog_soft_rf <- mse.cv.om_slog_soft_rf <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
  mpiw.cv.lcs_hard_rf <- cov.cv.lcs_hard_rf <- mse.cv.lcs_hard_rf <- array(NA,c(folds,nClusts, nCovars))
  mpiw.cv.lcs_soft_rf <- cov.cv.lcs_soft_rf <- mse.cv.lcs_soft_rf <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
  mpiw.cv.demos <- cov.cv.demos  <- mse.cv.demos <- array(NA,c(folds, nCovars - 1))


  # For adding in Sequence Metric Principal Components
  mpiw.cv.om_slog_hard_rf_3pc <- cov.cv.om_slog_hard_rf_3pc <- mse.cv.om_slog_hard_rf_3pc  <- array(NA,c(folds,nClusts, nCovars + 3)) 
  mpiw.cv.lcs_soft_rf_3pc  <- cov.cv.lcs_soft_rf_3pc  <- mse.cv.lcs_soft_rf_3pc <- array(NA,c(folds, nClusts, nCovars + nClusts + 3 - 2))
  mpiw.cv.om_trate_soft_rf_3pc  <- cov.cv.om_trate_soft_rf_3pc  <- mse.cv.om_trate_soft_rf_3pc <- array(NA,c(folds, nClusts, nCovars + nClusts + 3 - 2))
  mpiw.cv.om_slog_soft_rf_3pc  <- cov.cv.om_slog_soft_rf_3pc  <- mse.cv.om_slog_soft_rf_3pc <- array(NA,c(folds, nClusts, nCovars + nClusts + 3 - 2))
  mpiw.cv.harm_rf_3pc <- cov.cv.harm_rf_3pc <- mse.cv.harm_rf_3pc <- array(NA,c(folds,nHarms, nCovars + nHarms + 3 - 1 ))
  mpiw.cv.windows_rf_3pc <- cov.cv.windows_rf_3pc <- mse.cv.windows_rf_3pc <- array(NA,c(folds,nWindows, nCovars + nWindows + 3 - 1))


  # for tracking convergence / condition number (do not need the mtry dimension)
  con.cv.lcs_soft <- con.cv.om_trate_soft <- con.cv.om_slog_soft <- array(NA, c(folds, nClusts))

  # Data readying 
  rmetrics <- create_rmetrics(mvad.seq=mvad.seq)
  mvad_long <- create_long_data(mvad)
  basis <- create.bspline.basis(c(0, 35), nbasis = 6, norder = 4)
  year_state_counts <- create_year_state_counts(mvad)

  ### Demos
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx

    # Make principal components from sequence metrics pcs 
    seqmets_pcs <- create_seqmets_pcs(rmetrics, train_idx, test_idx, nSeqMetsPCs)
    train_seq_scores <- seqmets_pcs$train_pcs
    test_seq_scores <- seqmets_pcs$test_pcs
    
    # data frame of mvad demos split into training and testing 
    demos_data <- mvad_covars |> mutate(y = num_month_em_last_year)

    for (j in 1:(nCovars - 1)) {
      
      # Fit with different mtry 
      demos_output <- fit_rf(demos_data[train_idx,], demos_data[test_idx,], mtry = j)

      mse.cv.demos[i, j] <- demos_output$mse
      cov.cv.demos[i, j] <- demos_output$coverage
      mpiw.cv.demos[i, j] <- demos_output$mpiw
    }
  
    # CFDA 
    harm_data <- create_cfda_harms(mvad_long, train_idx, test_idx, basis=basis)
    train_harm <- harm_data$train_data
    test_harm <- harm_data$test_data

    for (j in 1:nHarms) {

      train_j_harm <- train_harm[, 1:c(12+j)]
      test_j_harm <- test_harm[, 1:c(12+j)]
      

      for (k in 1:(nCovars + j - 1)) {

        harm_fit <- fit_rf(train_j_harm,test_j_harm, mtry = k)
        mse.cv.harm_rf[i, j, k] <- harm_fit$mse
        cov.cv.harm_rf[i, j, k] <- harm_fit$coverage
        mpiw.cv.harm_rf[i, j, k] <- harm_fit$mpiw

        # Adding in Sequence Metrics 
        harm_fit_seq <- fit_rf(cbind(train_j_harm, train_seq_scores),cbind(test_j_harm, test_seq_scores), mtry = k)
        mse.cv.harm_rf_3pc[i, j, k] <- harm_fit_seq$mse
        cov.cv.harm_rf_3pc[i, j, k] <- harm_fit_seq$coverage
        mpiw.cv.harm_rf_3pc[i, j, k] <- harm_fit_seq$mpiw

        }
      }

    count_data <- create_counts_pcs(year_state_counts, mvad_covars, train_idx, test_idx, num_month_em_last_year)
    train_count_data <- count_data$train_data
    test_count_data <- count_data$test_data
    
    for (j in 1:nWindows) {
      train_j_wind <- train_count_data[,1:(12+j)]
      test_j_wind <- test_count_data[,1:(12+j)]

      for (k in 1:(nCovars + j - 1)) {
        wind_fit <- fit_rf(
          train_data=train_j_wind, test_data = test_j_wind, mtry=k)
        mse.cv.windows_rf[i,j,k] <- wind_fit$mse
        cov.cv.windows_rf[i, j, k] <- wind_fit$coverage
        mpiw.cv.windows_rf[i, j, k] <- wind_fit$mpiw
        
        # Add in Sequence Metrics 
        wind_fit_seq <- fit_rf(
          train_data=cbind(train_j_wind, train_seq_scores), test_data = cbind(test_j_wind, test_seq_scores), mtry=k)
        mse.cv.windows_rf_3pc[i,j,k] <- wind_fit_seq$mse
        cov.cv.windows_rf_3pc[i, j, k] <- wind_fit_seq$coverage
        mpiw.cv.windows_rf_3pc[i, j, k] <- wind_fit_seq$mpiw

      }
    }


    # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
    clusterward_hard <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")


    # hard coded for OM-Trate
    for (j in 2:nClusts) {

      # OM INDEL SLOG
      hard_cluster_data <- hard_cluster(
        clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, train_idx=train_idx, 
        test_idx=test_idx, dist_matrix = dists[[3]])

      train_om_hard <- hard_cluster_data$train_data
      test_om_hard <- hard_cluster_data$test_data
    
      # OM-Slog
      soft_cluster_data_om_slog <- soft_cluster(dist_matrix = dists[[3]], train_idx=train_idx, test_idx=test_idx, nClusts=j, 
        covars=mvad_covars, y = num_month_em_last_year, max_iters = max_soft_iters, fuzziness = fuzz_soft, 
         knn_soft_assign=knn_soft_assign, tol=tolerance)
      con.cv.om_slog_soft[i,j] <- soft_cluster_data_om_slog$converged

      train_om_slog_soft <- soft_cluster_data_om_slog$train_data
      test_om_slog_soft <- soft_cluster_data_om_slog$test_data

      
      # metrics data 
      train_om_hard_3pc <- cbind(train_om_hard, train_seq_scores)
      test_om_hard_3pc <- cbind(test_om_hard, test_seq_scores)

      train_om_slog_soft_3pc <- cbind(train_om_slog_soft, train_seq_scores)
      test_om_slog_soft_3pc <- cbind(test_om_slog_soft, test_seq_scores)
      

      for (k in 1:(nCovars + j - 2)) {

        if (k < (nCovars + 1)) {
          hard_slog_fit <- fit_rf(train_data=train_om_hard, test_data=test_om_hard, mtry=k)
          mse.cv.om_slog_hard_rf[i,j,k] <- hard_slog_fit$mse
          cov.cv.om_slog_hard_rf[i, j, k] <- hard_slog_fit$coverage
          mpiw.cv.om_slog_hard_rf[i, j, k] <- hard_slog_fit$mpiw
        
          hard_slog_3pc_fit <- fit_rf(train_data=train_om_hard_3pc, test_data=test_om_hard_3pc, mtry=k)
          mse.cv.om_slog_hard_rf_3pc[i,j,k] <- hard_slog_3pc_fit$mse
          cov.cv.om_slog_hard_rf_3pc[i,j,k] <- hard_slog_3pc_fit$coverage
          mpiw.cv.om_slog_hard_rf_3pc[i,j,k] <- hard_slog_3pc_fit$mpiw

        }


        # OM-Slog
        if (soft_cluster_data_om_slog$converged) {
          soft_om_slog_fit <- fit_rf(train_data=train_om_slog_soft, test_data=test_om_slog_soft, mtry=k)
          mse.cv.om_slog_soft_rf[i,j,k] <- soft_om_slog_fit$mse
          cov.cv.om_slog_soft_rf[i, j, k] <- soft_om_slog_fit$coverage
          mpiw.cv.om_slog_soft_rf[i, j, k] <- soft_om_slog_fit$mpiw

          soft_om_slog_fit_3pc <- fit_rf(train_data=train_om_slog_soft_3pc, test_data=test_om_slog_soft_3pc, mtry=k)
          mse.cv.om_slog_soft_rf_3pc[i,j,k] <- soft_om_slog_fit_3pc$mse
          cov.cv.om_slog_soft_rf_3pc[i,j,k] <- soft_om_slog_fit_3pc$coverage
          mpiw.cv.om_slog_soft_rf_3pc[i,j,k] <- soft_om_slog_fit_3pc$mpiw

        }

      }
    }
  }



  mse.comp <- list()
  mse.comp$om_trate_hard <- mse.cv.om_trate_hard_rf  
  mse.comp$om_trate_soft <- mse.cv.om_trate_soft_rf 
  mse.comp$om_slog_hard  <- mse.cv.om_slog_hard_rf 
  mse.comp$om_slog_soft  <- mse.cv.om_slog_soft_rf 
  mse.comp$lcs_hard <- mse.cv.lcs_hard_rf 
  mse.comp$lcs_soft <- mse.cv.lcs_soft_rf 
  mse.comp$windows <- mse.cv.windows_rf 
  mse.comp$harm <- mse.cv.harm_rf 
  mse.comp$demos <- mse.cv.demos
  mse.comp$om_trate_soft_3pc <-  mse.cv.om_trate_soft_rf_3pc 
  mse.comp$om_slog_soft_3pc <-  mse.cv.om_slog_soft_rf_3pc 
  mse.comp$lcs_soft_3pc <-  mse.cv.lcs_soft_rf_3pc 
  mse.comp$om_slog_hard_3pc <- mse.cv.om_slog_hard_rf_3pc 
  mse.comp$windows_3pc <- mse.cv.windows_rf_3pc
  mse.comp$harm_3pc <- mse.cv.harm_rf_3pc
  
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
  cov.comp$om_trate_soft_3pc <-  cov.cv.om_trate_soft_rf_3pc 
  cov.comp$om_slog_soft_3pc <-  cov.cv.om_slog_soft_rf_3pc 
  cov.comp$lcs_soft_3pc <-  cov.cv.lcs_soft_rf_3pc 
  cov.comp$om_slog_hard_3pc <- cov.cv.om_slog_hard_rf_3pc 
  cov.comp$windows_3pc <- cov.cv.windows_rf_3pc
  cov.comp$harm_3pc <- cov.cv.harm_rf_3pc
      
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
  mpiw.comp$om_trate_soft_3pc <- mpiw.cv.om_trate_soft_rf_3pc 
  mpiw.comp$om_slog_soft_3pc <- mpiw.cv.om_slog_soft_rf_3pc 
  mpiw.comp$lcs_soft_3pc <- mpiw.cv.lcs_soft_rf_3pc 
  mpiw.comp$om_slog_hard_3pc <- mpiw.cv.om_slog_hard_rf_3pc  
  mpiw.comp$windows_3pc <- mpiw.cv.windows_rf_3pc
  mpiw.comp$harm_3pc <- mpiw.cv.harm_rf_3pc

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



saveRDS(cv_comp, file="CV_NonLinearComp_Factor_w3pc_r_checks_more_soft.rds")

  
# Stop the cluster
stopCluster(cl)

  
  
  
  
  
  
  
  
  
  
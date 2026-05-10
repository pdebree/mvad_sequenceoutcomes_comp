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

nCVs <- 2

comp <- list()
comp$om_trate_hard <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$om_trate_soft <- array(NA,c(folds,nSoftClusts,nCovars + nSoftClusts - 2))
comp$om_slog_hard <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$om_slog_soft <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
comp$lcs_hard <- array(NA,c(folds,nClusts, nCovars + nClusts - 2))
comp$lcs_soft <- array(NA,c(folds,nSoftClusts, nCovars + nSoftClusts - 2))
comp$windows <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))
comp$harm <- array(NA,c(folds,nHarms, nCovars + nHarms - 1 ))

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
      for (k in 1:(nCovars + j - 1)) {
          # create training set based on number of PCs to include
        train_harm <- mvad_covars[train_idx, ] %>% 
          add_column(as_tibble(pcs.train[, 1:j, drop = FALSE])) %>% 
          mutate(y = num_month_em_last_year[train_idx])
        test_harm <- mvad_covars[test_idx, ] %>% 
          add_column(as_tibble(pcs.test[, 1:j, drop = FALSE])) %>% 
          mutate(y = num_month_em_last_year[test_idx])

        harm_fit <- fit_rf(train_harm, test_harm, mtry = k)
        mse.cv.harm_rf[i, j, k] <- harm_fit$mse
        cov.cv.harm_rf[i, j, k] <- harm_fit$cov
        mpiw.cv.harm_rf[i, j, k] <- harm_fit$mpiw


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
    
    colnames(train_mvad_windows)[1] <- "y"
    colnames(test_mvad_windows)[1] <- "y"
    
    for (j in 1:nWindows) {
      for (k in 1:(nCovars + j - 1)) {
        wind_fit <- fit_rf(
          train_data=train_mvad_windows[,1:(11+j)], test_data = test_mvad_windows, mtry=k)
        mse.cv.windows_rf[i,j,k] <- wind_fit$mse
        cov.cv.windows_rf[i, j, k] <- wind_fit$cov
        mpiw.cv.windows_rf[i, j, k] <- wind_fit$mpiw
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
      print(j)
      # Create hard clusters
      hard_cluster_data <- hard_cluster_onehot(
        clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, train_idx=train_idx, 
        test_idx=test_idx, dist_matrix = dists[[1]])

      train_om_hard <- hard_cluster_data$train_data
      test_om_hard <- hard_cluster_data$test_data
    
      # Create soft clusters
      soft_cluster_data <- soft_cluster(dist_matrix = dists[[1]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, y = num_month_em_last_year)

      train_om_soft <- soft_cluster_data$train_data
      test_om_soft <- soft_cluster_data$test_data

      # have to minus 2 because of indexing j from 2 (instead of 1 - can't look at 1 soft cluster)
      for (k in 1:(nCovars + j - 2)) {
        hard_om_fit <- fit_rf(train_data=train_om_hard, test_data=test_om_hard, mtry=k)
        mse.cv.om_trate_hard_rf[i,j,k] <- hard_om_fit$mse
        cov.cv.om_trate_hard_rf[i, j, k] <- hard_om_fit$cov
        mpiw.cv.om_trate_hard_rf[i, j, k] <- hard_om_fit$mpiw

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
      
      # Create hard clusters
      hard_cluster_data <- hard_cluster_onehot(
        clusterward = clusterward_hard, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, train_idx=train_idx, 
        test_idx=test_idx, dist_matrix = dists[[3]])

      train_om_hard <- hard_cluster_data$train_data
      test_om_hard <- hard_cluster_data$test_data
    
      # Create soft clusters
      soft_cluster_data <- soft_cluster(dist_matrix = dists[[3]], 
        train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year)

      train_om_soft <- soft_cluster_data$train_data
      test_om_soft <- soft_cluster_data$test_data

      for (k in 1:(nCovars + j - 2)) {
        hard_slog_fit <- fit_rf(train_data=train_om_hard, test_data=test_om_hard, mtry=k)
        mse.cv.om_slog_hard_rf[i,j,k] <- hard_slog_fit$mse
        cov.cv.om_slog_hard_rf[i, j, k] <- hard_slog_fit$cov
        mpiw.cv.om_slog_hard_rf[i, j, k] <- hard_slog_fit$mpiw
        
        # make sure we don't look at 25 columns (but we still want to look at mtry up to the number of soft clusters)
        # because it has multiple columns (one for each cluster)
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

    # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
    clusterward_hard <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")

    for (j in 2:nClusts) {
      
      cluster_data <- hard_cluster_onehot(clusterward = clusterward_hard, nClusts = j, covars=mvad_covars, y = num_month_em_last_year, train_idx = train_idx, test_idx = test_idx, dist_matrix = dists[[2]])
      
      train_lcs_hard <- cluster_data$train_data
      test_lcs_hard <- cluster_data$test_data

      soft_cluster_data <- soft_cluster(dist_matrix = dists[[2]], train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, y = num_month_em_last_year)

      train_lcs_soft <- soft_cluster_data$train_data
      test_lcs_soft <- soft_cluster_data$test_data
      
      for (k in 1:(nCovars + j - 2)) {

        hard_lcs_fit <- fit_rf(train_data=train_lcs_hard, test_data=test_lcs_hard, mtry=k)
        mse.cv.lcs_hard_rf[i,j,k] <- hard_lcs_fit$mse
        cov.cv.lcs_hard_rf[i, j, k] <- hard_lcs_fit$cov
        mpiw.cv.lcs_hard_rf[i, j, k] <- hard_lcs_fit$mpiw

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

  
  cov.comp <- list()
  cov.comp$om_trate_hard <- cov.cv.om_trate_hard_rf 
  cov.comp$om_trate_soft <- cov.cv.om_trate_soft_rf 
  cov.comp$om_slog_hard  <- cov.cv.om_slog_hard_rf 
  cov.comp$om_slog_soft  <- cov.cv.om_slog_soft_rf 
  cov.comp$lcs_hard <- cov.cv.lcs_hard_rf 
  cov.comp$lcs_soft <- cov.cv.lcs_soft_rf 
  cov.comp$windows <- cov.cv.windows_rf 
  cov.comp$harm <- cov.cv.harm_rf

  mpiw.comp <- list()
  mpiw.comp$om_trate_hard <-  mpiw.cv.om_trate_hard_rf 
  mpiw.comp$om_trate_soft <-  mpiw.cv.om_trate_soft_rf 
  mpiw.comp$om_slog_hard  <- mpiw.cv.om_slog_hard_rf 
  mpiw.comp$om_slog_soft  <-  mpiw.cv.om_slog_soft_rf 
  mpiw.comp$lcs_hard <- mpiw.cv.lcs_hard_rf 
  mpiw.comp$lcs_soft <- mpiw.cv.lcs_soft_rf 
  mpiw.comp$windows <- mpiw.cv.windows_rf 
  mpiw.comp$harm <-  mpiw.cv.harm_rf 

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
    
}

for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_NonLinearComp.rds")

  
# Stop the cluster
stopCluster(cl)

  
  
  
  
  
  
  
  
  
  
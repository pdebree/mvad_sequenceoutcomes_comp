# CFDA Linear Competition 
# First written for Multi-Level Models practicum in Fall 2025

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
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

  ### Cross Validation for Clustering Methods with OM-Transition Rate 
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx
    
    # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM-Trate
    clusterward_hard <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")
    
    # hard coded for OM-Trate
    for (j in 2:nClusts) {

      hard_cluster_data <- hard_cluster(
        clusterward=clusterward_hard, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, train_idx=train_idx, 
        test_idx=test_idx, dist_matrix = dists[[1]])

    
      hard_output <- fit_linear(hard_cluster_data$train_data, hard_cluster_data$test_data)



      # hard clustering
      mse.cv.om_trate_hard[i,j] <- hard_output$mse
      cov.cv.om_trate_hard[i,j] <- hard_output$coverage
      mpiw.cv.om_trate_hard[i,j] <- hard_output$mpiw
      
      #soft clustering - fanny 
      if (j < nSoftClusts) {
        soft_cluster_data <- soft_cluster(dist_matrix = dists[[1]], 
          train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
          y = num_month_em_last_year)

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
    
    # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM-Slog
    clusterward_hard <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")
    
    # hard coded for OM-Trate
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
      
      #soft clustering - fanny 
      if (j < nSoftClusts) {
        soft_cluster_data <- soft_cluster(dist_matrix = dists[[3]], 
          train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
          y = num_month_em_last_year)

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
    
    # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR LCS-DSS
    clusterward_hard <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")
    
    for (j in 2:nClusts) {

      hard_cluster_data <- hard_cluster(
        clusterward=clusterward_hard, nClusts=j, covars=mvad_covars, 
        y = num_month_em_last_year, train_idx=train_idx, 
        test_idx=test_idx, dist_matrix = dists[[2]])

      
      fit.lin <- lm(y~., data=hard_cluster_data$train_data) 
      preds.lin <- predict(fit.lin, newdata = hard_cluster_data$test_data, interval = "prediction", level = 0.95)
      mse.lin <- calc_mse(preds.lin[,"fit"], y_dat=hard_cluster_data$test_data[,"y"])
      

      hard_output <- fit_linear(hard_cluster_data$train_data, hard_cluster_data$test_data)

      # hard clustering
      mse.cv.lcs_hard[i,j] <- hard_output$mse
      cov.cv.lcs_hard[i,j] <- hard_output$coverage
      mpiw.cv.lcs_hard[i,j] <- hard_output$mpiw
      
      #soft clustering - fanny 
      if (j < nSoftClusts) {
        soft_cluster_data <- soft_cluster(dist_matrix = dists[[2]], 
          train_idx=train_idx, test_idx=test_idx, nClusts=j, covars=mvad_covars, 
          y = num_month_em_last_year)

        soft_output <- fit_linear(soft_cluster_data$train_data, soft_cluster_data$test_data)
        
        if (soft_cluster_data$converged) {
          mse.cv.lcs_soft[i,j] <- soft_output$mse
          cov.cv.lcs_soft[i,j] <- soft_output$coverage
          mpiw.cv.lcs_soft[i,j] <- soft_output$mpiw
        }
      }
    }
  }

  ### Year-State Count Windows (PCA)

  # Make our windows of data - first three years
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
      
    pca_windows_train <- prcomp(x=mvad_windows[12:27][train_idx,], center=TRUE, scale=TRUE)
    train_scores <- pca_windows_train$x
    test_scores <- predict(pca_windows_train, newdata = mvad_windows[12:27][test_idx, ]) 

    colnames(train_scores) <- paste0("PC", 1:ncol(train_scores))
    colnames(test_scores) <- paste0("PC", 1:ncol(test_scores))
    
    # add principal components to demographic data (no longer need original windows)
    train_mvad_windows <- cbind(mvad_windows[1:11][train_idx,], as.data.frame(train_scores))
    test_mvad_windows <- cbind(mvad_windows[1:11][test_idx,], as.data.frame(test_scores))
    
    for (j in 1:nWindows) {

      train_j_wind <- train_mvad_windows[, 1:c(11+j)] |> mutate(y = num_month_em_last_year[train_idx])
      test_j_wind <- test_mvad_windows[, 1:c(11+j)] |> mutate(y = num_month_em_last_year[test_idx])

      wind_fit <- fit_linear(train_j_wind, test_j_wind)
      mse.cv.windows[i,j] <- wind_fit$mse
      cov.cv.windows[i,j] <- wind_fit$coverage
      mpiw.cv.windows[i,j] <- wind_fit$mpiw
    }
  }




  ### CFDA 

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
    pcs.train <- fmca.train$pc
    nComps <- ncol(pcs.train)
    
    colnames(pcs.train) <- paste0("PC",1:nComps)
    pcs.test <- predict(fmca.train,newdata=mvad_long_test,method="parallel",nCores=7)
    colnames(pcs.test) <- paste0("PC",1:nComps)
    
    #it's relatively easy to run this with 5 or 10 harmonics - or in the equation, below
    train_harm <- mvad_covars[train_idx,] %>% add_column(as_tibble(pcs.train)) %>% mutate(y=num_month_em_last_year[train_idx])
    test_harm <- mvad_covars[test_idx,] %>% add_column(as_tibble(pcs.test)) %>% mutate(y=num_month_em_last_year[test_idx])
    
    for (j in 1:nHarms) {

      train_harm <- mvad_covars[train_idx,] %>% add_column(as_tibble(pcs.train)[,1:j]) %>% mutate(y=num_month_em_last_year[train_idx])
      test_harm <- mvad_covars[test_idx,] %>% add_column(as_tibble(pcs.test)[,1:j]) %>% mutate(y=num_month_em_last_year[test_idx])
      
      harm_fit <- fit_linear(train_harm, test_harm)
      
      mse.cv.harm[i,j] <- harm_fit$mse
      cov.cv.harm[i,j] <- harm_fit$coverage
      mpiw.cv.harm[i,j] <- harm_fit$mpiw

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

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
}

# Populate global matrices from return list 



for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_LinearComp.rds")

  
# Stop the cluster
stopCluster(cl)



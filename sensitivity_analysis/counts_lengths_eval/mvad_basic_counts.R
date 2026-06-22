# Basic Counts - Linear and NonLinear Runs

# Pippi de Bree
# gcd2056@nyu.edu
# 2026-05-28


library(tidyverse)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)

source("seqout_utils.R")


folds <- 5
nWindows <- 6 

nCVs <- 20

comp <- list()
comp$bow_windows <- array(NA,c(folds,nWindows))

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

  nWindows <- 6 


  ### Arrays to hold RMSEs from different runs 
  mse.cv.windows <- cov.cv.windows <- mpiw.cv.windows <- array(NA,c(folds,nWindows))

  ### Year-State Count Windows (PCA)

  # Make our windows of data - first three years
  mvad_states <- mvad[,c(1,15:50)]

  # makes year and id combinations with counts for each state
  id_states <- mvad_states %>% pivot_longer(cols=colnames(mvad_states[2:37]), names_to="output") %>% count(id, value)

  # make full table with 16 columns (one for each year/id combo) + 1 for id.
  mvad_states_wide <- id_states %>% 
    pivot_wider(id_cols=c(id), 
                names_from=value, values_from=n, values_fill = 0) %>% dplyr::select(-id)

  mvad_windows <- cbind(mvad_covars, mvad_states_wide)


  for (i in 1:folds) {

    cat("Fold Number ",i,"\n")
    test_idx <- idx == i  
    train_idx <- !test_idx
      
    pca_windows_train <- prcomp(x=mvad_windows[12:17][train_idx,], center=TRUE, scale=TRUE)
    train_scores <- pca_windows_train$x
    test_scores <- predict(pca_windows_train, newdata = mvad_windows[12:17][test_idx, ]) 

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


  mse.comp <- list()
  mse.comp$windows <- mse.cv.windows

  cov.comp <- list()
  cov.comp$windows <- cov.cv.windows

  mpiw.comp <- list()
  mpiw.comp$windows <- mpiw.cv.windows
  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
}

for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_LinearComp_basicwindows.rds")

  
# Stop the cluster
# stopCluster(cl)


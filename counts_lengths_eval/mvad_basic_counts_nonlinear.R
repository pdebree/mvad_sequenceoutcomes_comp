# 

# Pippi de Bree
# gcd2056@nyu.edu
# 2026-05-28


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
library(doRNG)

source("seqout_utils.R")


# Detect cores allocated by Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# Register the cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(seed = 123)


folds <- 5
nWindows <- 6 
nCovars <- 11 

nCVs <- 20

comp <- list()
comp$windows <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))

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


  nWindows <- 6

  # Arrays for holding outcomes fits 
  # soft clusters - 2 (because we never look at index=1 clusters and indexing is always from 1)
  mpiw.cv.windows_rf <- cov.cv.windows_rf <- mse.cv.windows_rf <- array(NA,c(folds,nWindows, nCovars + nWindows - 1))

  ### Windows 
  mvad_states <- mvad[,c(1,15:50)]

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

    # pull out training and testing for only the windows
    train_windows <- mvad_windows[12:17][train_idx, ]
    test_windows <- mvad_windows[12:17][test_idx, ]
    
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

  mse.comp <- list()
  mse.comp$windows <- mse.cv.windows_rf 
  
  cov.comp <- list()
  cov.comp$windows <- cov.cv.windows_rf 

  mpiw.comp <- list()
  mpiw.comp$windows <- mpiw.cv.windows_rf 

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
    
}

for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_NonLinearComp_basicwindows.rds")

  
  
    
# Stop the cluster
stopCluster(cl)

  
  
  
  
  
  
  
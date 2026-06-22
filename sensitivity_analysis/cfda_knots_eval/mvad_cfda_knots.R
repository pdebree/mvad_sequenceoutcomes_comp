# 

# Pippi de Bree
# gcd2056@nyu.edu
# 2026-05-29


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
nHarms <- 25
nCVs <- 20

comp <- list()
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
  nHarms <- 25


  ### Arrays to hold RMSEs from different runs 
  mse.cv.harm <- cov.cv.harm  <- mpiw.cv.harm <- array(NA,c(folds,nHarms))


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
  basis <- create.bspline.basis(c(0, M), nbasis = 10, norder = 4)


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
  mse.comp$harm <- mse.cv.harm


  cov.comp <- list()
  cov.comp$harm <- cov.cv.harm


  mpiw.comp <- list()
  mpiw.comp$harm <- mpiw.cv.harm

  list(mpiw.comp=mpiw.comp, cov.comp=cov.comp, mse.comp=mse.comp, cv_index=cv_index)
}

# Populate global matrices from return list 



for(result in results) {
  n <- result$cv_index 
  cv_comp$mse[[n]] <- result$mse.comp
  cv_comp$cov[[n]] <- result$cov.comp
  cv_comp$mpiw[[n]] <- result$mpiw.comp
}


saveRDS(cv_comp, file="CV_LinearComp_cfda10.rds")

  
# Stop the cluster
# stopCluster(cl)



# Simulated Data Evaluation of Life Course Methods

library(tidyverse)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)

source("seqout_utils.R")

# Detect cores allocated by Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# Register the cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)




# Make the ordering of the methods strictly as follows: 
# 1 - cfda
# 2 - windows 
# 3 - om-trate (hard)
# 4 - om-slog (hard)
# 5 - lcs (hard)
# 6 - om-trate (soft)
# 7 - om-slog (soft)
# 8 - lcs (soft)

nMethods <- 8
nComps <- 6
nSoft <- 4
nSets <- 150
fuzz_soft <- 1.75
max_soft_iters <- 2000
tolerance <- 1e-10
knn_soft_assign <- 15

folds <- 10
cv_idx <- rep(1:folds,each=floor(900/folds))
    
sim_types <- c(
  "dar_med_seqs", 
  "easy_seqs", 
  "hard_seqs", 
  "med_seqs", 
  "very_hard_seqs"
)

sim_type <- "hard_seqs"

# Data is stored in directory based on the difficulty types.
# The "hard_seqs" directory contains the following files)
file_names <- c(
  "simDat33_concord.rdata", 
  "simDat33_indep_semi.rdata",
  "simDat33_indep.rdata",
  "simDat33_semi_concord.rdata",
  "simDat33_semi.rdata")

# have to strictly follow the ordering of the methods 
# Overall performance (best in each method) for each 
# data set 
mse.sims <- list()
mse.sims$concord <- array(NA, c(nMethods, nSets))
mse.sims$indep_semi <- array(NA, c(nMethods, nSets))
mse.sims$indep <- array(NA, c(nMethods, nSets))
mse.sims$semi_concord <- array(NA, c(nMethods, nSets))
mse.sims$semi <- array(NA, c(nMethods, nSets))


cov.sims <- list()
cov.sims$concord <- array(NA, c(nMethods, nSets))
cov.sims$indep_semi <- array(NA, c(nMethods, nSets))
cov.sims$indep <- array(NA, c(nMethods, nSets))
cov.sims$semi_concord <- array(NA, c(nMethods, nSets))
cov.sims$semi <- array(NA, c(nMethods, nSets))


mpiw.sims <- list()
mpiw.sims$concord <- array(NA, c(nMethods, nSets))
mpiw.sims$indep_semi <- array(NA, c(nMethods, nSets))
mpiw.sims$indep <- array(NA, c(nMethods, nSets))
mpiw.sims$semi_concord <- array(NA, c(nMethods, nSets))
mpiw.sims$semi <- array(NA, c(nMethods, nSets))

# locations for training and testing sets - first 900 
# are training, second 900 are testing. 
set_train_idx <- c(rep(TRUE, 900), rep(FALSE, 900))
set_test_idx <- !set_train_idx


# data structure for holding the best number of 
# components by performance method. 
best_n_comps <- list()
best_n_comps$concord <- array(NA, c(nMethods, nSets))
best_n_comps$indep_semi <- array(NA, c(nMethods, nSets))
best_n_comps$indep <- array(NA, c(nMethods, nSets))
best_n_comps$semi_concord <- array(NA, c(nMethods, nSets))
best_n_comps$semi <- array(NA, c(nMethods, nSets))

set_conv_train_track <- list()
set_conv_train_track$om_trate <- set_conv_train_track$om_slog <- set_conv_train_track$lcs <- array(FALSE, c(folds, nSoft, 2))

train_conv_tracker <- list()
train_conv_tracker$concord <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$indep_semi <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$indep <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$semi_concord <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$semi <- replicate(nSets, set_conv_train_track, simplify = FALSE)


set_conv_test_track <- list()
set_conv_test_track$om_trate <- set_conv_test_track$om_slog <- set_conv_test_track$lcs <- array(NA, 2)

test_conv_tracker <- list()
test_conv_tracker$concord <- replicate(nSets, set_conv_test_track, simplify = FALSE)
test_conv_tracker$indep_semi  <- replicate(nSets, set_conv_test_track, simplify = FALSE)
test_conv_tracker$indep <- replicate(nSets, set_conv_test_track, simplify = FALSE)
test_conv_tracker$semi_concord <- replicate(nSets, set_conv_test_track, simplify = FALSE)
test_conv_tracker$semi <- replicate(nSets, set_conv_test_track, simplify = FALSE)


# Flattened version (potentially faster)
task_grid <- expand.grid(
  file_name = file_names,
  set_idx = 1:nSets,
  stringsAsFactors = FALSE
)

results <- foreach(m = 1:nrow(task_grid), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {
  
  # load current data 
  load(paste0("sim_data/", sim_type, "/", task_grid$file_name[m]))

  if (!exists("df_long_list")) {
        stop(paste("df_long_list not found in file:", task_grid$file_name[m]))
    }


  n <- task_grid$set_idx[m]
  sim_group <- sub("^[^_]*_(.*)\\..*$", "\\1", task_grid$file_name[m])

  cur_df_long <- df_long_list[[n]]
    
  # pick set to do cross validation on 
  data_wide <- cur_df_long %>% group_by(ID) %>% 
    mutate(y = y[time=="T1"]) %>% ungroup() %>%  # make every y value equal to the "first" value  
    pivot_wider(id_cols=c(ID, Class, y), names_from=time, values_from=state) %>% 
    dplyr::select(-Class, -ID) # drop Class variable and ID is the same as the index 
  
  # create distance matrices (first column is Y value) 
  dists <- create_dists(data.seq=seqdef(data_wide, 2:13, 
                    states=c("A", "B", "C")))
  
  # split data into train and testing for the set
  set_train_data <- data_wide[set_train_idx,]
  set_test_data <- data_wide[set_test_idx,]
  

  # Re-initialize for every cv (every set)
  mse.cv <- array(NA, c(folds, nComps, nMethods))
  
  # CFDA training and testing data
  cfda_data_long <- cur_df_long %>% mutate(
    time= as.numeric(str_sub(time, start = 2, end = -1)),
    state=as.factor(state)) %>%
    rename(id=ID) %>%
    dplyr::select(-Class, -y)

  cfda_long_set_train <- cfda_data_long %>% filter(id %in% 1:900)
  cfda_long_set_test <- cfda_data_long %>% filter(id %in% 901:1800)


  # Windows training and testing data
  set_data_windows <- cfda_data_long %>% mutate(
    window = if_else(time < 7, 1, 2)) %>%
    count(id, window, state) %>%
    pivot_wider(
      id_cols=c(id), names_from=c(state,window), values_from=n, values_fill = 0)

  
  # soft clustering convergence tracking 
  set_conv_train_track <- list()
  set_conv_train_track$om_trate <- set_conv_train_track$om_slog <- set_conv_train_track$lcs <- array(FALSE, c(folds, nSoft, 2))

  # initialize basis functions - can do this generally
  basis <- create.bspline.basis(c(1, 12), nbasis = 3, norder = 3)
  
  # cross validation
  for(i in 1:folds) {
      cv_test_idx <- cv_idx == i
      cv_train_idx <- !cv_test_idx
      
      # get fold train and test sets, from training set.
      y.cv_train <- set_train_data[cv_train_idx, "y"]
      y.cv_test <- set_train_data[cv_test_idx, "y"]
      
    
      # Hard clusterings - will "cut" depending on number of clusters. 
      clusterward_trate_hard <- agnes(
        dists[[1]][set_train_idx, set_train_idx][cv_train_idx,cv_train_idx],
          diss=TRUE, method="ward")
      clusterward_slog_hard <- agnes(
        dists[[3]][set_train_idx, set_train_idx][cv_train_idx,cv_train_idx],
          diss=TRUE, method="ward")
      clusterward_lcs_hard <- agnes(
        dists[[2]][set_train_idx, set_train_idx][cv_train_idx,cv_train_idx],
          diss=TRUE, method="ward")


      # Creating CFDA harmonics 
      fmca.train <- compute_optimal_encoding(cfda_long_set_train %>% filter(id %in% seq(1:900)[cv_train_idx]),
                                            basis,verbose=F)
      pcs.train <- fmca.train$pc
        
      nCFDAComps <- ncol(pcs.train)

      
      pcs.test <- predict(fmca.train,newdata=
                            cfda_long_set_train %>% filter(id %in% seq(1:900)[cv_test_idx]),verbose=F)
    
      cv_train_harm <- data.frame(y = y.cv_train) %>% add_column(as_tibble(pcs.train) %>% 
        set_names(paste0("PC", 1:ncol(pcs.train))))
      cv_test_harm <- data.frame(y = y.cv_test) %>% add_column(as_tibble(pcs.test) %>%
        set_names(paste0("PC", 1:ncol(pcs.test))))
      colnames(cv_test_harm) <- colnames(cv_train_harm) # ensure names match

      # 
      # Windows Counts 
      pca_windows_train <- prcomp(x=set_data_windows[2:7][set_train_idx, ][cv_train_idx,])

      # add principal components to demographic data (no longer need original windows)
      cv_train_windows <- cbind(as.data.frame(y.cv_train), pca_windows_train$x)
      cv_test_windows <- cbind(as.data.frame(y.cv_test),
                              predict(
                                pca_windows_train,
                                newdata = set_data_windows[2:7][set_train_idx, ][cv_test_idx,]) )

    
      for(j in 1:nComps) {
        
        if (j > 1) {
          # ***************
          # hard clustering. 
          om_trate_clusters <- hard_cluster(clusterward = clusterward_trate_hard, nClusts = j,  y = set_train_data[["y"]], 
            train_idx = cv_train_idx, test_idx = cv_test_idx, dist_matrix = dists[[1]][set_train_idx, set_train_idx])

          om_slog_clusters <- hard_cluster(clusterward = clusterward_slog_hard, nClusts = j,  y = set_train_data[["y"]], 
            train_idx = cv_train_idx, test_idx = cv_test_idx,dist_matrix = dists[[3]][set_train_idx, set_train_idx])
          
          lcs_clusters <- hard_cluster(clusterward = clusterward_lcs_hard, nClusts = j, y = set_train_data[["y"]], 
            train_idx = cv_train_idx, test_idx = cv_test_idx, dist_matrix = dists[[2]][set_train_idx, set_train_idx])


          # hard is method "1" - fit model and calculate mse
          mse.cv[i, j, 3] <- fit_linear(
            train_data = om_trate_clusters$train_data,
            test_data = om_trate_clusters$test_data)$mse
          mse.cv[i, j, 4] <- fit_linear(
            train_data = om_slog_clusters$train_data,
            test_data = om_slog_clusters$test_data)$mse
          mse.cv[i, j, 5] <- fit_linear(
            train_data = lcs_clusters$train_data,
            test_data = lcs_clusters$test_data)$mse
          
        # **********************
        # soft clustering - we have to consider potential convergence problems
          if (j <= nSoft) {
            om_trate_soft_train <- soft_cluster_sim(dists[[1]][set_train_idx, set_train_idx],
              cv_train_idx, cv_test_idx, nClusts=j, fuzziness = fuzz_soft, y=set_train_data[["y"]], max_iters = max_soft_iters, 
              knn_soft_assign=knn_soft_assign, tol=tolerance
            )

            set_conv_train_track[["om_trate"]][i,j,] <- c(om_trate_soft_train$converged[[1]], om_trate_soft_train$safe_condition_no)

            if (om_trate_soft_train$converged && om_trate_soft_train$safe_condition_no) {
              # hard is method "2" - fit model and calculate mse 
              mse.cv[i, j, 6] <- fit_linear(om_trate_soft_train$train_data, om_trate_soft_train$test_data)$mse
            } 

            om_slog_soft_train <- soft_cluster_sim(dists[[3]][set_train_idx, set_train_idx],
              cv_train_idx, cv_test_idx, nClusts=j, fuzziness = fuzz_soft, y=set_train_data[["y"]], max_iters = max_soft_iters,
              knn_soft_assign=knn_soft_assign, tol=tolerance
            )

            set_conv_train_track[["om_slog"]][i,j,] <- c(om_slog_soft_train$converged[[1]], om_slog_soft_train$condition_no)


            if (om_slog_soft_train$converged && om_slog_soft_train$condition_no) {
              # hard is method "2" - fit model and calculate mse 
              mse.cv[i, j, 7] <- fit_linear(om_slog_soft_train$train_data, om_slog_soft_train$test_data)$mse
            } 

            lcs_soft_train <- soft_cluster_sim(dists[[2]][set_train_idx, set_train_idx],
              cv_train_idx, cv_test_idx, nClusts=j, fuzziness = fuzz_soft, y=set_train_data[["y"]], max_iters = max_soft_iters, 
              knn_soft_assign=knn_soft_assign, tol=tolerance
            )

            set_conv_train_track[["lcs"]][i,j,] <- c(lcs_soft_train$converged[[1]], lcs_soft_train$condition_no)

            if (lcs_soft_train$converged && lcs_soft_train$safe_condition_no) {
              # hard is method "2" - fit model and calculate mse 
              mse.cv[i, j, 8] <- fit_linear(lcs_soft_train$train_data, lcs_soft_train$test_data)$mse
            } 
          }
        } 
      
        # ****************
        # CFDA 
        
        # y value is the first column, so we index to j + 1 to get j elements for this round.
        # j = 2 means we look at the first PC (so clusters we get 2 but PCs we get 1 for j=2)
        # "3" is CFDA
        mse.cv[i, j, 1] <- fit_linear(cv_train_harm[,1:(j+1)], cv_test_harm[,1:(j+1)])$mse

        # # ****************
        # # Windows

        # # can't do more than 6 PCs because there are only 6 dimensions
        # # (should really do less than 6)
        if (j < 7) {
          mse.cv[i, j, 2] <- fit_linear(cv_train_windows[,1:(j+1)], cv_test_windows[,1:(j+1)])$mse
        }
      }
  }
  

  mse_cfda <- mse_wind <- mse_om_trate_hard <- NA
  mse_om_slog_hard <- mse_lcs_hard <- mse_om_trate_soft <- NA
  mse_om_slog_soft <- mse_lcs_soft <- NA

  cov_cfda <- cov_wind <- cov_om_trate_hard <- NA
  cov_om_slog_hard <- cov_lcs_hard <- cov_om_trate_soft <- NA
  cov_om_slog_soft <- cov_lcs_soft <- NA

  mpiw_cfda <- mpiw_wind <- mpiw_om_trate_hard <- NA
  mpiw_om_slog_hard <- mpiw_lcs_hard <- mpiw_om_trate_soft <- NA
  mpiw_om_slog_soft <- mpiw_lcs_soft <- NA

  
  # Fitting by best parameters
  best_cfda_harms <- which.min(apply(mse.cv[,,1], 2, mean, na.rm = TRUE))
  best_windows_pcs <- which.min(apply(mse.cv[,,2], 2, mean, na.rm = TRUE))
  best_trate_hard_clusts <- which.min(apply(mse.cv[,,3], 2, mean, na.rm = TRUE))
  best_slog_hard_clusts <- which.min(apply(mse.cv[,,4], 2, mean, na.rm = TRUE))
  best_lcs_hard_clusts <- which.min(apply(mse.cv[,,5], 2, mean, na.rm = TRUE))
  
  # If we get non-convergence (which corresponds to an NA value for that entry in 
  # mse.cv) we ignore it. We base this on the assumption that non-convergence is not happening
  # more than twice in a group of folds for a number of soft clusters
  best_trate_soft_clusts <- safe_which_min(apply(mse.cv[,,6], 2, mean, na.rm = TRUE))
  best_slog_soft_clusts <- safe_which_min(apply(mse.cv[,,7], 2, mean, na.rm = TRUE))
  best_lcs_soft_clusts <- safe_which_min(apply(mse.cv[,,8], 2, mean, na.rm = TRUE))
  
  # Refitting models based on best training performance (number of components)
  y.train_set <- set_train_data[, "y"]
  y.test_set <- set_test_data[, "y"]

  clusterward_om_trate_hard_set <- agnes(
    dists[[1]][set_train_idx, set_train_idx],
    diss=TRUE, method="ward")
  clusterward_om_slog_hard_set <- agnes(
    dists[[3]][set_train_idx, set_train_idx],
    diss=TRUE, method="ward")
  clusterward_lcs_hard_set <- agnes(
    dists[[2]][set_train_idx, set_train_idx],
    diss=TRUE, method="ward")

  # Create clusterings based on best hard performance
  om_trate_clusters_set <- hard_cluster(clusterward = clusterward_om_trate_hard_set, 
    nClusts = best_trate_hard_clusts,  y = data_wide[["y"]], train_idx = set_train_idx, test_idx = set_test_idx, dist_matrix = dists[[1]])
  om_slog_clusters_set <- hard_cluster(clusterward = clusterward_om_slog_hard_set,
    nClusts = best_slog_hard_clusts, y = data_wide[["y"]], train_idx = set_train_idx, test_idx = set_test_idx, dist_matrix = dists[[3]])
  lcs_clusters_set <- hard_cluster(clusterward = clusterward_lcs_hard_set,
    nClusts = best_lcs_hard_clusts, y = data_wide[["y"]], train_idx = set_train_idx, test_idx = set_test_idx, dist_matrix = dists[[2]])

  hard_trate_set_fit <- fit_linear(
    train_data=om_trate_clusters_set$train_data, 
    test_data=om_trate_clusters_set$test_data)
  mse_om_trate_hard <- hard_trate_set_fit$mse
  cov_om_trate_hard <- hard_trate_set_fit$coverage
  mpiw_om_trate_hard <- hard_trate_set_fit$mpiw

  hard_slog_set_fit <- fit_linear(
    train_data=om_slog_clusters_set$train_data, 
    test_data=om_slog_clusters_set$test_data)
  mse_om_slog_hard <- hard_slog_set_fit$mse
  cov_om_slog_hard <- hard_slog_set_fit$coverage
  mpiw_om_slog_hard <- hard_slog_set_fit$mpiw

  hard_lcs_set_fit <- fit_linear(
    train_data=lcs_clusters_set$train_data, 
    test_data=lcs_clusters_set$test_data)
  mse_lcs_hard <- hard_lcs_set_fit$mse
  cov_lcs_hard <- hard_lcs_set_fit$coverage
  mpiw_lcs_hard <- hard_lcs_set_fit$mpiw

  
  # om trate soft 
  # if there is no training convergence, mse is NA 
  if (!is.na(best_trate_soft_clusts)) {
    om_trate_soft_set <- soft_cluster_sim(dists[[1]],
            set_train_idx, set_test_idx, nClusts=best_trate_soft_clusts, fuzziness = fuzz_soft, y=data_wide[["y"]], max_iters = max_soft_iters, 
          knn_soft_assign=knn_soft_assign, tol=tolerance)
    
    om_trate_soft_conv <- c(TRUE, om_trate_soft_set$converged[[1]] == 1, om_trate_soft_set$condition_no)
    
    if (om_trate_soft_set$converged && om_trate_soft_set$safe_condition_no) {
      soft_trate_set_fit <- fit_linear(train_data=om_trate_soft_set$train_data, test_data=om_trate_soft_set$test_data)
      mse_om_trate_soft <- soft_trate_set_fit$mse
      cov_om_trate_soft <- soft_trate_set_fit$coverage
      mpiw_om_trate_soft <- soft_trate_set_fit$mpiw
    } else {
      mse_om_trate_soft <- cov_om_trate_soft <- mpiw_om_trate_soft <- NA
    }

  } else {
    mse_om_trate_soft <- cov_om_trate_soft <- mpiw_om_trate_soft <- NA
    om_trate_soft_conv <- c(0, 0, 0)
  }

  # om-slog soft 
  if (!is.na(best_slog_soft_clusts)) {
    om_slog_soft_set <- soft_cluster_sim(dists[[3]],
            set_train_idx, set_test_idx, nClusts=best_slog_soft_clusts, fuzziness = fuzz_soft, y=data_wide[["y"]], max_iters = max_soft_iters, 
          knn_soft_assign=knn_soft_assign, tol=tolerance)
    
    om_slog_soft_conv <- c(TRUE, om_slog_soft_set$converged[[1]], om_slog_soft_set$condition_no)
    
    if (om_slog_soft_set$converged && om_slog_soft_set$safe_condition_no) {
      soft_slog_set_fit <- fit_linear(train_data=om_slog_soft_set$train_data, test_data=om_slog_soft_set$test_data)
      mse_om_slog_soft <- soft_slog_set_fit$mse
      cov_om_slog_soft <- soft_slog_set_fit$coverage
      mpiw_om_slog_soft <- soft_slog_set_fit$mpiw
    } else {
      mse_om_slog_soft <- cov_om_slog_soft <- mpiw_om_slog_soft <- NA
    }

  } else {
    mse_om_slog_soft <- cov_om_slog_soft <- mpiw_om_slog_soft <- NA
    om_slog_soft_conv <- c(0, 0, 0)
  }


  # lcs soft
  if (!is.na(best_lcs_soft_clusts)) {
    lcs_soft_set <- soft_cluster_sim(dists[[2]],
            set_train_idx, set_test_idx, nClusts=best_lcs_soft_clusts, fuzziness = fuzz_soft, y=data_wide[["y"]], max_iters = max_soft_iters, 
          knn_soft_assign=knn_soft_assign, tol=tolerance)
    
    lcs_soft_conv <- c(TRUE, lcs_soft_set$converged[[1]], lcs_soft_set$condition_no)
    
    if (lcs_soft_set$converged && lcs_soft_set$safe_condition_no) {
      soft_lcs_set_fit <- fit_linear(train_data=lcs_soft_set$train_data, test_data=lcs_soft_set$test_data)
      mse_lcs_soft <- soft_lcs_set_fit$mse
      cov_lcs_soft <- soft_lcs_set_fit$coverage
      mpiw_lcs_soft <- soft_lcs_set_fit$mpiw
    } else {
      mse_lcs_soft <- cov_lcs_soft <- mpiw_lcs_soft <- NA
    }
  } else {
    mse_lcs_soft <- cov_lcs_soft <- mpiw_lcs_soft <- NA
    lcs_soft_conv <- c(0, 0, 0)
  }


  # CFDA - refit
  fmca.train_set <- compute_optimal_encoding(cfda_long_set_train, basis,verbose=F)
  pcs.train_set <- fmca.train_set$pc
  pcs.test_set <- predict(fmca.train_set,newdata=cfda_long_set_test,verbose=F)

  set_train_harm <- data.frame(y = y.train_set) %>% add_column(as_tibble(pcs.train_set[,1:best_cfda_harms, drop = FALSE]))
  set_test_harm <- data.frame(y = y.test_set) %>% add_column(as_tibble(pcs.test_set[,1:best_cfda_harms, drop = FALSE]))
  colnames(set_test_harm) <- colnames(set_train_harm) # ensure names match

  cfda_set_fit <- fit_linear(set_train_harm, set_test_harm)
  mse_cfda <- cfda_set_fit$mse
  cov_cfda <- cfda_set_fit$coverage
  mpiw_cfda <- cfda_set_fit$mpiw
  

  # Windows - refit
  pca_windows_train_set <- prcomp(x=set_data_windows %>% filter(id %in% 1:900) %>% dplyr::select(-id))
  set_train_scores <- pca_windows_train_set$x
  set_test_scores <- predict(pca_windows_train_set,
                            newdata = set_data_windows %>% filter(id %in% 901:1800) %>% dplyr::select(-id))
  
  
  wind_fit <- fit_linear(
    train_data=cbind(as.data.frame(y.train_set), set_train_scores[, 1:best_windows_pcs, drop=FALSE]),
    test_data=cbind(as.data.frame(y.test_set), set_test_scores[, 1:best_windows_pcs, drop=FALSE]))   
  
  mse_wind <- wind_fit$mse
  cov_wind <- wind_fit$coverage
  mpiw_wind <- wind_fit$mpiw
  
# return section (for this one file-set combination)
result <- list(
      sim_group = sim_group,
      set_idx = n,
      mses = c(mse_cfda, mse_wind, mse_om_trate_hard, mse_om_slog_hard, mse_lcs_hard, mse_om_trate_soft, 
        mse_om_slog_soft, mse_lcs_soft),
      covs = c(cov_cfda, cov_wind, cov_om_trate_hard, cov_om_slog_hard, cov_lcs_hard, cov_om_trate_soft, 
        cov_om_slog_soft, cov_lcs_soft),
      mpiw = c(mpiw_cfda, mpiw_wind, mpiw_om_trate_hard, mpiw_om_slog_hard, mpiw_lcs_hard, mpiw_om_trate_soft, 
        mpiw_om_slog_soft, mpiw_lcs_soft),
      best_k = c(best_cfda_harms, best_windows_pcs, best_trate_hard_clusts, best_slog_hard_clusts,
        best_lcs_hard_clusts, best_trate_soft_clusts,best_slog_soft_clusts, best_lcs_soft_clusts), 
      set_conv_train_track = set_conv_train_track,
      set_conv_test_track = list(
        om_trate = om_trate_soft_conv, 
        om_slog =  om_slog_soft_conv,
        lcs = lcs_soft_conv)
    )
}

# Populate global matrices from return list 
for(result in results) {
  g <- result$sim_group
  n <- result$set_idx
  
  mse.sims[[g]][, n] <- result$mses
  cov.sims[[g]][,n] <- result$covs
  mpiw.sims[[g]][, n] <- result$mpiw
  best_n_comps[[g]][, n] <- result$best_k
  train_conv_tracker[[g]][[n]][["om_trate"]] <- result$set_conv_train_track$om_trate
  train_conv_tracker[[g]][[n]][["om_slog"]] <- result$set_conv_train_track$om_slog
  train_conv_tracker[[g]][[n]][["lcs"]] <- result$set_conv_train_track$lcs

  test_conv_tracker[[g]][[n]] <- result$set_conv_test_track

}


full_results <- list()
full_results$mse <- mse.sims
full_results$cov <- cov.sims
full_results$mpiw <- mpiw.sims
full_results$best_n_comps <- best_n_comps
full_results$train_conv_tracker <- train_conv_tracker
full_results$test_conv_tracker <- test_conv_tracker

# Save data from full run 
saveRDS(full_results, paste0(sim_type,"_results.rds"))

stopCluster(cl)

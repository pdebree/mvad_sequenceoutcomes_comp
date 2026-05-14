# CFDA Linear Competition 
# First written for Multi-Level Models practicum in Fall 2025

library(tidyverse)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)
library(ranger)
library(tuneRanger)
library(gt)
source("seqout_utils.R")


# Detect cores allocated by Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# Register the cluster
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 123) # for reproducability
registerDoParallel(cl)

nCVs <- 20
folds <- 5
nSeqPcs <- 12
nCovars <- 11
nClusts <- 25
nSoftClusts <- 13


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

task_vec <- data.frame(cv_index = 1:nCVs)

mse.cv.rmets_rf <- array(NA,c(folds, nSeqPcs, nCovars + nSeqPcs - 1))
cv_seqs <- replicate(nCVs, mse.cv.rmets_rf , simplify = FALSE)

seq_results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras", "ranger", "tuneRanger", "gt")) %dopar% {

  cv_index <- task_vec$cv_index[m]

  ### Training set up 
  folds <- 5
  nrecs <- nrow(mvad[1])

  # Set up training and test split so it is the same across the competition 
  idx.orig <- rep(1:folds,each=floor(nrecs/folds))
  if (nrecs %% folds != 0) idx.orig <- c(idx.orig,1:(nrecs %% folds))
  idx <- sample(idx.orig) #shuffle

  # Creates a random shuffle of the indices to be used in the the train-test split and folds. 
  shuffle_index <- sample(1:nrecs, nrecs)
  ids <- sort(unique(mvad$id))

  mse.cv.rmets_rf <- array(NA,c(folds, nSeqPcs, nCovars + nSeqPcs - 1))

  ## Clustering + Sequence Metrics 

  mvad_states <- mvad[,15:50]

  # encodings for the statebadness 
  st_alphabet <- alphabet(mvad.seq) 
  # Example: If alphabet is "A", "B", "C"
  st_prec_values <- c(1, -1, -1, 2, -1, -1) # A=1 (low badness), C=3 (high badness)
  names(st_prec_values) <- st_alphabet

  mvad_rmetrics <- mvad_states %>% mutate(
    spells=seqindic(mvad.seq, "dlgth")$Dlgth, 
    visited_states=seqindic(mvad.seq, "visited")$Visited, 
    num_of_trans = seqindic(mvad.seq,"trans")$Trans, 
    mean_spell_dur = seqindic(mvad.seq,"meand")$MeanD, 
    # pedantic - this one pulled out a "seqivardur" "numeric" datatype (not sure 
    # how it did both, so I have to force it to be numeric)
    sd_spell_dur = as.numeric(seqindic(mvad.seq,"dustd")$Dustd), 
    # Diversity I
    entropy = seqindic(mvad.seq,"entr")$Entr, 
    # more interested in states than spells (see paper)
    dss_subs = seqindic(mvad.seq,"nsubs")$Nsubs, 
    complexity = seqindic(mvad.seq,"cplx")$Cplx, 
    # could look at other turbulence measures 
    turbulence = seqindic(mvad.seq,"turb")$Turb, 
    badness = seqibad(seqdata = mvad.seq,stprec = st_prec_values), 
    degradation = seqidegrad(seqdata = mvad.seq, stprec = st_prec_values), 
    insecurity = seqinsecurity(seqdata = mvad.seq, stprec = st_prec_values))


  # Training set up 
  rmetrics <- mvad_rmetrics[,37:48] 
  #nSeqPcs <- ncol(rmetrics)

  # go through folds in repeat
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx
    
    train_seqmets <- rmetrics[train_idx, ]
    test_seqmets <- rmetrics[test_idx, ]
    
    pca_comps_train <- prcomp(x=train_seqmets, center=TRUE, scale=TRUE)
    
    train_scores <- pca_comps_train$x
    test_scores <- predict(pca_comps_train, newdata = test_seqmets) 
    
    for (j in 1:nSeqPcs) {
      train_rmets <- cbind(num_month_em_last_year[train_idx], mvad_covars[train_idx,], train_scores[,1:j,drop = FALSE])
      test_rmets <- cbind(num_month_em_last_year[test_idx], mvad_covars[test_idx,], test_scores[,1:j,drop = FALSE])
    
      colnames(train_rmets)[1] <- "y"
      colnames(test_rmets)[1] <- "y"
      
      for (k in 1:(nCovars + j - 1)) {
        mets_fit <- fit_rf(train_data=train_rmets, test_data=test_rmets, mtry=k)
        mse.cv.rmets_rf[i,j,k] <- mets_fit$mse
      }
    }
  }
  list(cv_index = cv_index, mses = mse.cv.rmets_rf)
}

for(seq_result in seq_results) {
  n <- seq_result$cv_index
  cv_seqs[[n]] <- seq_result$mses
}


seqs_list <- simplify2array(cv_seqs)
# do mean and then min because we want to find mean performance by mtry and then find the 
# number of components that has the best values for all mtry combos
mean_mse_vec <- apply(apply(seqs_list, c(2, 3), mean, na.rm = TRUE), 1, min, na.rm=TRUE)

best_n_seq_mets_pc <- safe_which_min(mean_mse_vec)

mse.seq_clusts <- list()
mse.seq_clusts$om_trate_hard <- array(NA,c(folds,nClusts,nCovars + best_n_seq_mets_pc + nClusts - 2))
mse.seq_clusts$om_trate_soft <-  array(NA,c(folds,nClusts,nCovars + best_n_seq_mets_pc  + nSoftClusts - 2))
mse.seq_clusts$om_slog_hard <- array(NA,c(folds,nClusts,nCovars  + best_n_seq_mets_pc + nClusts - 2))
mse.seq_clusts$om_slog_soft <- array(NA,c(folds,nClusts,nCovars  + best_n_seq_mets_pc + nSoftClusts - 2))
mse.seq_clusts$lcs_hard <- array(NA,c(folds,nClusts,nCovars + best_n_seq_mets_pc  + nClusts - 2))
mse.seq_clusts$lcs_soft <- array(NA,c(folds,nClusts,nCovars  + best_n_seq_mets_pc + nSoftClusts - 2))

cv_seq_clusts <- replicate(nCVs, mse.seq_clusts, simplify = FALSE)



clust_results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {

  cv_index <- task_vec$cv_index[m]
  
  # Number of each method to check
  nClusts <- 25 
  nSoftClusts <- 13 
  nSeqPcs <- 12
  fuzz_soft <- 1.5

  mse.seq_clusts <- list()
  mse.seq_clusts$om_trate_hard <- array(NA,c(folds,nClusts,nCovars + best_n_seq_mets_pc + nClusts - 2))
  mse.seq_clusts$om_trate_soft <-  array(NA,c(folds,nClusts,nCovars + best_n_seq_mets_pc  + nSoftClusts - 2))
  mse.seq_clusts$om_slog_hard <- array(NA,c(folds,nClusts,nCovars  + best_n_seq_mets_pc + nClusts - 2))
  mse.seq_clusts$om_slog_soft <- array(NA,c(folds,nClusts,nCovars  + best_n_seq_mets_pc + nSoftClusts - 2))
  mse.seq_clusts$lcs_hard <- array(NA,c(folds,nClusts,nCovars + best_n_seq_mets_pc  + nClusts - 2))
  mse.seq_clusts$lcs_soft <- array(NA,c(folds,nClusts,nCovars  + best_n_seq_mets_pc + nSoftClusts - 2))


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


  ### Sequence PCs + Clustering Methods 
  # need to find the best performance with just the sequence PCs then add in the 

  mvad_states <- mvad[,15:50]

  # encodings for the statebadness 
  st_alphabet <- alphabet(mvad.seq) 
  # Example: If alphabet is "A", "B", "C"
  st_prec_values <- c(1, -1, -1, 2, -1, -1) # A=1 (low badness), C=3 (high badness)
  names(st_prec_values) <- st_alphabet

  mvad_rmetrics <- mvad_states %>% mutate(
    spells=seqindic(mvad.seq, "dlgth")$Dlgth, 
    visited_states=seqindic(mvad.seq, "visited")$Visited, 
    num_of_trans = seqindic(mvad.seq,"trans")$Trans, 
    mean_spell_dur = seqindic(mvad.seq,"meand")$MeanD, 
    # pedantic - this one pulled out a "seqivardur" "numeric" datatype (not sure 
    # how it did both, so I have to force it to be numeric)
    sd_spell_dur = as.numeric(seqindic(mvad.seq,"dustd")$Dustd), 
    # Diversity I
    entropy = seqindic(mvad.seq,"entr")$Entr, 
    # more interested in states than spells (see paper)
    dss_subs = seqindic(mvad.seq,"nsubs")$Nsubs, 
    complexity = seqindic(mvad.seq,"cplx")$Cplx, 
    # could look at other turbulence measures 
    turbulence = seqindic(mvad.seq,"turb")$Turb, 
    badness = seqibad(seqdata = mvad.seq,stprec = st_prec_values), 
    degradation = seqidegrad(seqdata = mvad.seq, stprec = st_prec_values), 
    insecurity = seqinsecurity(seqdata = mvad.seq, stprec = st_prec_values))

  rmetrics <- mvad_rmetrics[,37:48] 
  nSeqPcs <- ncol(rmetrics) 

  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx
    
    # using rmetrics from above
    train_rmetrics <- rmetrics[train_idx, ]
    test_rmetrics <- rmetrics[test_idx, ]
    
    pca_comps_train <- prcomp(x=train_rmetrics, center=TRUE, scale=TRUE)
    
    train_scores <- pca_comps_train$x[,1:best_n_seq_mets_pc]
    test_scores <- predict(pca_comps_train, newdata = test_rmetrics)[,1:best_n_seq_mets_pc]
        
    # Do 6 different combinations of clustering and sequence PCS
    # create hierarchical clustering
    clusterward_hard_trate <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")
    clusterward_hard_lcs <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")
    clusterward_hard_slog <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")

    # create clusterings 


    for (j in 2:nClusts) {

      om_trate_clustered_hard <- hard_cluster_onehot(clusterward = clusterward_hard_trate, nClusts = j, 
        covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[1]])
      om_trate_hard_train <- cbind(om_trate_clustered_hard$train_data, train_scores) 
      om_trate_hard_test <- cbind(om_trate_clustered_hard$test_data, test_scores) 
      
      # om-slog hard
    om_slog_clustered_hard <- hard_cluster_onehot(clusterward = clusterward_hard_slog, nClusts = j, 
        covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[3]])
    om_slog_hard_train <- cbind(om_slog_clustered_hard$train_data, train_scores) 
    om_slog_hard_test <- cbind(om_slog_clustered_hard$test_data, test_scores)
      
    # lcs hard
    lcs_clustered_hard <- hard_cluster_onehot(clusterward = clusterward_hard_lcs, nClusts = j, 
        covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[2]])
    lcs_hard_train <- cbind(lcs_clustered_hard$train_data, train_scores)
    lcs_hard_test <- cbind(lcs_clustered_hard$test_data, test_scores) 
      
    # om-trate soft
      if (j < nSoftClusts) {
        om_trate_clustered_soft <- soft_cluster(dists[[1]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        om_trate_soft_train <- cbind(om_trate_clustered_soft$train_data, train_scores)
        om_trate_soft_test <- cbind(om_trate_clustered_soft$test_data, test_scores) 

        # om-slog soft
        om_slog_clustered_soft <- soft_cluster(dists[[3]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        om_slog_soft_train <- cbind(om_slog_clustered_soft$train_data, train_scores)
        om_slog_soft_test <- cbind(om_slog_clustered_soft$test_data, test_scores)

        # lcs soft 
        lcs_clustered_soft <- soft_cluster(dists[[2]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        lcs_soft_train <- cbind(lcs_clustered_soft$train_data, train_scores) 
        lcs_soft_test <- cbind(lcs_clustered_soft$test_data, test_scores) 
        
      }

      for (k in 1:(nCovars + best_n_seq_mets_pc + j - 2)) {

        hard_seq_trate_fit <- fit_rf(
          train_data = om_trate_hard_train, 
          test_data = om_trate_hard_test, mtry=k)
        mse.seq_clusts[["om_trate_hard"]][i, j,k] <- hard_seq_trate_fit$mse
        
        hard_seq_slog_fit <- fit_rf(
          train_data = om_slog_hard_train, 
          test_data = om_slog_hard_test, mtry=k)
        mse.seq_clusts[["om_slog_hard"]][i, j,k] <- hard_seq_slog_fit$mse

        hard_seq_lcs_fit <-  fit_rf(
          train_data = lcs_hard_train, 
          test_data = lcs_hard_test, mtry=k)
        mse.seq_clusts[["lcs_hard"]][i, j,k] <- hard_seq_lcs_fit$mse

        
        if (j < nSoftClusts) {
        
          if (om_trate_clustered_soft$converged) {
            soft_seq_trate_fit <- fit_rf(
              train_data = om_trate_soft_train, 
              test_data = om_trate_soft_test, mtry=k)
            mse.seq_clusts[["om_trate_soft"]][i,j,k] <- soft_seq_trate_fit$mse
          }
        
          if (om_slog_clustered_soft$converged) {
            soft_seq_slog_fit <- fit_rf(
              train_data = om_slog_soft_train, 
              test_data = om_slog_soft_test, mtry=k)
            mse.seq_clusts[["om_slog_soft"]][i,j,k] <- soft_seq_slog_fit$mse

          }

          if (lcs_clustered_soft$converged) {
            soft_seq_lcs_fit <- fit_rf(
              train_data = lcs_soft_train, 
              test_data = lcs_soft_test, mtry=k)
            mse.seq_clusts[["lcs_soft"]][i,j,k] <- soft_seq_lcs_fit$mse
          }
        }
      }
    }
  }

  list(
    cv_index = cv_index, 
    cv_mses = mse.seq_clusts
  )
}

for(clust_result in clust_results) {
  n <- clust_result$cv_index
  cv_seq_clusts[[n]] <- clust_result$cv_mses
}



get_method_mses_mtry_frame <- function(seq_mses, method) {

  all_mses <- lapply(seq_mses, function(cv_data) {
    cv_data[[method]]
    })
  
  combined_4d <- simplify2array(all_mses)

  # mean across folds and cvs
  mean_mses <- sqrt(apply(combined_4d, c(2, 3), mean, na.rm = TRUE)) |> 
    as.data.frame() |> rownames_to_column(var = "Index") 
  colnames(mean_mses)[-1] <- 1:(ncol(mean_mses) - 1)

  df <- mean_mses |> pivot_longer(cols = -Index, names_to = "mtry", values_to = "rmse")
  df$method <- method

  return(df)
}

seq_clusts_mses <- get_method_mses_mtry_frame(cv_seq_clusts, "om_trate_hard")



for (method in  c("om_trate_soft", "om_slog_hard", "om_slog_soft", "lcs_hard", "lcs_soft")) {

  seq_clusts_mses <- rbind(seq_clusts_mses, get_method_mses_mtry_frame(cv_seq_clusts, method))
}


saveRDS(seq_clusts_mses, file="CV_SeqMetsNonLinear.rds")
 
# Stop the cluster
stopCluster(cl)
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

nCVs <- 2
folds <- 10 
nSeqPcs <- 12

task_vec <- data.frame(cv_index = 1:nCVs)

mse.mets  <- array(NA,c(folds,nSeqPcs))
cv_seqs <- replicate(nCVs, mse.mets, simplify = FALSE)

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


seq_results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {

  cv_index <- task_vec$cv_index[m]


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
  nSoftClusts <- 13 
  nSeqPcs <- 12
  fuzz_soft <- 1.5

  mse.cv.mets  <- array(NA,c(folds,nSeqPcs))
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

  # Loop over folds to find best performance 
  for (i in 1:folds) {
    cat("Fold Number ",i,"\n")
    test_idx <- idx == i
    train_idx <- !test_idx

    
    pca_comps_train <- prcomp(x=rmetrics[train_idx, ], center=TRUE, scale=TRUE)
    
    train_scores <- pca_comps_train$x
    test_scores <- predict(pca_comps_train, newdata =rmetrics[test_idx, ]) 
    
    train_mvad_rmets <- cbind(mvad_covars[train_idx,], as.data.frame(train_scores))
    test_mvad_rmets <- cbind(mvad_covars[test_idx,], as.data.frame(test_scores))
    
    # find best performance in fold 
    for (j in 1:nSeqPcs) {

        # This is fairly hard coded to the form of mvad_covars0. 
      train_seq <- train_mvad_rmets[,1:(11+j)] %>% mutate(y=num_month_em_last_year[train_idx])
      test_seq <- test_mvad_rmets[,1:(11+j)] %>% mutate(y=num_month_em_last_year[test_idx])

      seqs_fit <- fit_linear(train_seq, test_seq)

      mse.cv.mets[i,j] <- seqs_fit$mse
    }
  }

  list(cv_index=cv_index, mse = mse.cv.mets)
}


for(seq_result in seq_results) {
  n <- seq_result$cv_index
  cv_seqs[[n]] <- seq_result$mse

}


seq_3d <- simplify2array(cv_seqs)
seq_rmses <- sqrt(apply(seq_3d, 2, mean))
best_n_seq_mets_pc <- which.min(seq_rmses)

  
mse.seq_clusts <- list()
mse.seq_clusts$om_trate_hard <- array(NA,c(folds,nClusts))
mse.seq_clusts$om_trate_soft <-  array(NA,c(folds,nClusts))
mse.seq_clusts$om_slog_hard <- array(NA,c(folds,nClusts))
mse.seq_clusts$om_slog_soft <- array(NA,c(folds,nClusts))
mse.seq_clusts$lcs_hard <- array(NA,c(folds,nClusts))
mse.seq_clusts$lcs_soft <- array(NA,c(folds,nClusts))

cv_seq_clusts <- replicate(nCVs, mse.seq_clusts, simplify = FALSE)



clust_results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {

  cv_index <- task_vec$cv_index[m]
  
  # Old way for handling best way to do it 
  # data frame to hold best combinations - long makes more sense 
  mse.seq_clusts <- list()
  mse.seq_clusts$om_trate_hard <- array(NA,c(folds,nClusts))
  mse.seq_clusts$om_trate_soft <-  array(NA,c(folds,nClusts))
  mse.seq_clusts$om_slog_hard <- array(NA,c(folds,nClusts))
  mse.seq_clusts$om_slog_soft <- array(NA,c(folds,nClusts))
  mse.seq_clusts$lcs_hard <- array(NA,c(folds,nClusts))
  mse.seq_clusts$lcs_soft <- array(NA,c(folds,nClusts))


    # Number of each method to check
  nClusts <- 25 
  nSoftClusts <- 13 
  nSeqPcs <- 12
  fuzz_soft <- 1.5

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
  nSoftClusts <- 13 
  nSeqPcs <- 12
  fuzz_soft <- 1.5


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
    
    train_pcs <- as.data.frame(pca_comps_train$x[,1:best_n_seq_mets_pc])
    test_pcs <- as.data.frame(predict(pca_comps_train, newdata = test_rmetrics)[,1:best_n_seq_mets_pc])

    colnames(train_pcs) <- paste0("SeqPC", 1:best_n_seq_mets_pc)
    colnames(test_pcs)  <- paste0("SeqPC", 1:best_n_seq_mets_pc)
        
    # Do 6 different combinations of clustering and sequence PCS
    # create hierarchical clustering
    clusterward_hard_trate <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")
    clusterward_hard_lcs <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")
    clusterward_hard_slog <- agnes(dists[[3]][train_idx,train_idx], diss=TRUE, method="ward")

    # create clusterings 

    for (j in 2:nClusts) {

      om_trate_clustered_hard <- hard_cluster(clusterward = clusterward_hard_trate, nClusts = j, 
        covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[1]])
      om_trate_hard_train <- cbind(om_trate_clustered_hard$train_data, train_pcs)
      om_trate_hard_test <- cbind(om_trate_clustered_hard$test_data, test_pcs)
      
      mse.seq_clusts[["om_trate_hard"]][i, j] <- fit_linear(
          train_data = om_trate_hard_train, 
          test_data = om_trate_hard_test)$mse
      
      # om-slog hard
      om_slog_clustered_hard <- hard_cluster(clusterward = clusterward_hard_slog, nClusts = j, 
          covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[3]])
      
      om_slog_hard_train <- cbind(om_slog_clustered_hard$train_data, train_pcs)
      om_slog_hard_test <- cbind(om_slog_clustered_hard$test_data, test_pcs)
        
      mse.seq_clusts[["om_slog_hard"]][i, j] <- fit_linear(
            train_data = om_slog_hard_train, 
            test_data = om_slog_hard_test)$mse

        
      # lcs hard
      lcs_clustered_hard <- hard_cluster(clusterward = clusterward_hard_lcs, nClusts = j, 
          covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[2]])
      lcs_hard_train <- cbind(lcs_clustered_hard$train_data, train_pcs)
      lcs_hard_test <- cbind(lcs_clustered_hard$test_data, test_pcs)
        

      mse.seq_clusts[["lcs_hard"]][i, j] <- fit_linear(train_data = lcs_hard_train, test_data = lcs_hard_test)$mse
        
    # om-trate soft
      if (j < nSoftClusts) {
        om_trate_clustered_soft <- soft_cluster(dists[[1]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        
        if (om_trate_clustered_soft$converged) {
          om_trate_soft_train <- cbind(om_trate_clustered_soft$train_data, train_pcs)
          om_trate_soft_test <- cbind(om_trate_clustered_soft$test_data, test_pcs)
          mse.seq_clusts[["om_trate_soft"]][i,j] <- fit_linear(
              train_data = om_trate_soft_train, 
              test_data = om_trate_soft_test)$mse
          }
        
        # om-slog soft
        om_slog_clustered_soft <- soft_cluster(dists[[3]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        if (om_slog_clustered_soft$converged) {
          om_slog_soft_train <- cbind(om_slog_clustered_soft$train_data, train_pcs)
          om_slog_soft_test <- cbind(om_slog_clustered_soft$test_data, test_pcs)
          mse.seq_clusts[["om_slog_soft"]][i,j] <- fit_linear(
              train_data = om_slog_soft_train, 
              test_data = om_slog_soft_test)$mse
        }

        # lcs soft 
        lcs_clustered_soft <- soft_cluster(dists[[2]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        if (lcs_clustered_soft$converged) {
          lcs_soft_train <- cbind(lcs_clustered_soft$train_data, train_pcs)
          lcs_soft_test <- cbind(lcs_clustered_soft$test_data, test_pcs)
          mse.seq_clusts[["lcs_soft"]][i,j] <- fit_linear(
              train_data = lcs_soft_train, 
              test_data = lcs_soft_test)$mse
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
  n <- seq_result$cv_index
  cv_seq_clusts[[n]] <- clust_result$cv_mses
}



get_method_mses <- function(seq_mses, method) {

  all_mses <- lapply(seq_mses, function(cv_data) {
    cv_data[[method]]
    })
  
  combined_4d <- simplify2array(all_mses)
  # mean across folds and cvs
  mean_mses <- sqrt(apply(combined_4d, 2, mean, na.rm = TRUE))
  return(mean_mses)
}


seq_clusts_mses <- data.frame(
  Index = 1:nClusts,
  n_seq_mets = best_n_seq_mets_pc,
  om_trate_hard = get_method_mses(cv_seq_clusts, "om_trate_hard"), 
  om_trate_soft = get_method_mses(cv_seq_clusts, "om_trate_soft"), 
  om_slog_hard = get_method_mses(cv_seq_clusts, "om_slog_hard"),
  om_slog_soft = get_method_mses(cv_seq_clusts, "om_slog_soft"),
  lcs_hard = get_method_mses(cv_seq_clusts, "lcs_hard"),
  lcs_soft = get_method_mses(cv_seq_clusts, "lcs_soft"))


saveRDS(seq_clusts_mses, file="CV_SeqMetsLinear.rds")
 
# Stop the cluster
stopCluster(cl)



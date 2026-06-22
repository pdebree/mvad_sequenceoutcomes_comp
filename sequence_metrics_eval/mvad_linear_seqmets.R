# CFDA Linear Competition 
# First written for Multi-Level Models practicum in Fall 2025

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)
library(doRNG)
source("seqout_utils.R")


# # Detect cores allocated by Slurm
# n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# # Register the cluster
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
# registerDoRNG(seed = 123)

nCVs <- 20
folds <- 5
nSeqPcs <- 12

# Number of each method to check
nClusts <- 25 
nSoftClusts <- 13 
fuzz_soft <- 1.5

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
avg_fold_mses <- apply(seq_3d, c(2,3), mean)
best_nseq_by_fold <- apply(avg_fold_mses, 2, safe_which_min)
best_n_seq_mets_pc <- safe_mode(best_nseq_by_fold)
avg_cv_mses <- apply(seq_3d, 2, mean)


comp <- list() 
comp$om_trate_hard <- array(NA,c(folds,nClusts))
comp$om_trate_soft <- array(NA,c(folds,nClusts))
comp$om_slog_hard <- array(NA,c(folds,nClusts))
comp$om_slog_soft <- array(NA,c(folds,nClusts))
comp$lcs_hard <- array(NA,c(folds,nClusts))
comp$lcs_soft <- array(NA,c(folds,nClusts))


cv_comp <- list()
cv_comp$mse <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$cov <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$mpiw <- replicate(nCVs, comp, simplify = FALSE)
cv_comp$n_seq_mets_modal <- best_n_seq_mets_pc
cv_comp$n_seq_mets_dist <- best_nseq_by_fold
cv_comp$seq_mses <- seq_3d


clust_results <- foreach(m = 1:nrow(task_vec), .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {

  cv_index <- task_vec$cv_index[m]
  
  mse.cv.om_trate_hard <- cov.cv.om_trate_hard <- mpiw.cv.om_trate_hard <- array(NA,c(folds,nClusts))
  mse.cv.om_trate_soft <- cov.cv.om_trate_soft <- mpiw.cv.om_trate_soft <- array(NA,c(folds,nClusts))
  mse.cv.om_slog_hard <- cov.cv.om_slog_hard <-  mpiw.cv.om_slog_hard <- array(NA,c(folds,nClusts))
  mse.cv.om_slog_soft <- cov.cv.om_slog_soft <- mpiw.cv.om_slog_soft <- array(NA,c(folds,nClusts))
  mse.cv.lcs_hard <- cov.cv.lcs_hard <- mpiw.cv.lcs_hard <- array(NA,c(folds,nClusts))
  mse.cv.lcs_soft <- cov.cv.lcs_soft <- mpiw.cv.lcs_soft <- array(NA,c(folds,nClusts))


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
      
      om_trate_hard_fit <- fit_linear(
          train_data = om_trate_hard_train, 
          test_data = om_trate_hard_test)
      
      mse.cv.om_trate_hard[i,j] <- om_trate_hard_fit$mse
      cov.cv.om_trate_hard[i,j] <- om_trate_hard_fit$coverage
      mpiw.cv.om_trate_hard[i,j] <- om_trate_hard_fit$mpiw
      

      # om-slog hard
      om_slog_clustered_hard <- hard_cluster(clusterward = clusterward_hard_slog, nClusts = j, 
          covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[3]])
      
      om_slog_hard_train <- cbind(om_slog_clustered_hard$train_data, train_pcs)
      om_slog_hard_test <- cbind(om_slog_clustered_hard$test_data, test_pcs)
        
      om_slog_hard_fit <- fit_linear(
            train_data = om_slog_hard_train, 
            test_data = om_slog_hard_test)

      mse.cv.om_slog_hard[i,j] <- om_slog_hard_fit$mse
      cov.cv.om_slog_hard[i,j] <- om_slog_hard_fit$coverage
      mpiw.cv.om_slog_hard[i,j] <- om_slog_hard_fit$mpiw
        
      # lcs hard
      lcs_clustered_hard <- hard_cluster(clusterward = clusterward_hard_lcs, nClusts = j, 
          covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[2]])
      lcs_hard_train <- cbind(lcs_clustered_hard$train_data, train_pcs)
      lcs_hard_test <- cbind(lcs_clustered_hard$test_data, test_pcs)
        

      lcs_hard_fit <- fit_linear(train_data = lcs_hard_train, test_data = lcs_hard_test)

      mse.cv.lcs_hard[i,j] <- lcs_hard_fit$mse
      cov.cv.lcs_hard[i,j] <- lcs_hard_fit$coverage
      mpiw.cv.lcs_hard[i,j] <- lcs_hard_fit$mpiw


    # om-trate soft
      if (j < nSoftClusts) {
        om_trate_clustered_soft <- soft_cluster(dists[[1]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        
        if (om_trate_clustered_soft$converged) {
          om_trate_soft_train <- cbind(om_trate_clustered_soft$train_data, train_pcs)
          om_trate_soft_test <- cbind(om_trate_clustered_soft$test_data, test_pcs)
          om_trate_soft_fit <- fit_linear(
              train_data = om_trate_soft_train, 
              test_data = om_trate_soft_test)
          mse.cv.om_trate_soft[i,j] <- om_trate_soft_fit$mse
          cov.cv.om_trate_soft[i,j] <- om_trate_soft_fit$coverage
          mpiw.cv.om_trate_soft[i,j] <- om_trate_soft_fit$mpiw
            
        }
        

        # om-slog soft
        om_slog_clustered_soft <- soft_cluster(dists[[3]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        if (om_slog_clustered_soft$converged) {
          om_slog_soft_train <- cbind(om_slog_clustered_soft$train_data, train_pcs)
          om_slog_soft_test <- cbind(om_slog_clustered_soft$test_data, test_pcs)
          om_slog_soft_fit <- fit_linear(
              train_data = om_slog_soft_train, 
              test_data = om_slog_soft_test)
          mse.cv.om_slog_soft[i,j] <- om_slog_soft_fit$mse
          cov.cv.om_slog_soft[i,j] <- om_slog_soft_fit$coverage
          mpiw.cv.om_slog_soft[i,j] <- om_slog_soft_fit$mpiw
          
        }

        # lcs soft 
        lcs_clustered_soft <- soft_cluster(dists[[2]], train_idx, test_idx = test_idx,
          nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
        if (lcs_clustered_soft$converged) {
          lcs_soft_train <- cbind(lcs_clustered_soft$train_data, train_pcs)
          lcs_soft_test <- cbind(lcs_clustered_soft$test_data, test_pcs)
          lcs_soft_fit <- fit_linear(
              train_data = lcs_soft_train, 
              test_data = lcs_soft_test)
          mse.cv.lcs_soft[i,j] <- lcs_soft_fit$mse
          cov.cv.lcs_soft[i,j] <- lcs_soft_fit$coverage
          mpiw.cv.lcs_soft[i,j] <- lcs_soft_fit$mpiw
        }
      }
    }
  }

    mse.comp <- list()
  mse.comp$om_trate_hard <- mse.cv.om_trate_hard 
  mse.comp$om_trate_soft <- mse.cv.om_trate_soft 
  mse.comp$om_slog_hard  <- mse.cv.om_slog_hard 
  mse.comp$om_slog_soft  <- mse.cv.om_slog_soft 
  mse.comp$lcs_hard <- mse.cv.lcs_hard
  mse.comp$lcs_soft <- mse.cv.lcs_soft


  cov.comp <- list()
  cov.comp$om_trate_hard <- cov.cv.om_trate_hard 
  cov.comp$om_trate_soft <- cov.cv.om_trate_soft 
  cov.comp$om_slog_hard  <- cov.cv.om_slog_hard 
  cov.comp$om_slog_soft  <- cov.cv.om_slog_soft 
  cov.comp$lcs_hard <- cov.cv.lcs_hard
  cov.comp$lcs_soft <- cov.cv.lcs_soft


  mpiw.comp <- list()
  mpiw.comp$om_trate_hard <- mpiw.cv.om_trate_hard 
  mpiw.comp$om_trate_soft <- mpiw.cv.om_trate_soft 
  mpiw.comp$om_slog_hard  <- mpiw.cv.om_slog_hard 
  mpiw.comp$om_slog_soft  <- mpiw.cv.om_slog_soft 
  mpiw.comp$lcs_hard <- mpiw.cv.lcs_hard
  mpiw.comp$lcs_soft <- mpiw.cv.lcs_soft


  list(
    cv_index = cv_index, 
    cv_mse = mse.comp, 
    cv_cov = cov.comp, 
    cv_mpiw = mpiw.comp
  )
}

for(clust_result in clust_results) {
  n <- clust_result$cv_index
  cv_comp$mse[[n]] <- clust_result$cv_mse
  cv_comp$cov[[n]] <- clust_result$cv_cov
  cv_comp$mpiw[[n]] <- clust_result$cv_mpiw

}

saveRDS(cv_comp, file="CV_SeqMetsLinear.rds")
 
# Stop the cluster
stopCluster(cl)



# CFDA Linear Competition 
# First written for Multi-Level Models practicum in Fall 2025

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
source("seqout_utils.R")

set.seed(11)

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
    cov.cv.mets[i,j] <- seqs_fit$coverage
    mpiw.cv.mets[i,j] <- seqs_fit$mpiw
  }
}

seq_mets_pc_means_lm <- sqrt(apply(mse.cv.mets, 2, mean))
best_n_seq_mets_pc <- which.min(seq_mets_pc_means_lm)

# Old way for handling best way to do it 
# data frame to hold best combinations - long makes more sense 
mse.seq_clusts <- list()
mse.seq_clusts$om_trate_hard <- array(NA,c(folds,nClusts))
mse.seq_clusts$om_trate_soft <-  array(NA,c(folds,nClusts))
mse.seq_clusts$om_slog_hard <- array(NA,c(folds,nClusts))
mse.seq_clusts$om_slog_soft <- array(NA,c(folds,nClusts))
mse.seq_clusts$lcs_hard <- array(NA,c(folds,nClusts))
mse.seq_clusts$lcs_soft <- array(NA,c(folds,nClusts))


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

  # om_trate_hard - 

  for (j in 2:nClusts) {

    om_trate_clustered_hard <- hard_cluster(clusterward = clusterward_hard_trate, nClusts = j, 
      covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[1]])
    om_trate_hard_train <- cbind(om_trate_clustered_hard$train_data, as.data.frame(train_scores))
    om_trate_hard_test <- cbind(om_trate_clustered_hard$test_data, as.data.frame(test_scores))
    
    mse.seq_clusts[["om_trate_hard"]][i, j] <- fit_linear(
        train_data = om_trate_hard_train, 
        test_data = om_trate_hard_test)$mse
    
    # om-slog hard
    om_slog_clustered_hard <- hard_cluster(clusterward = clusterward_hard_slog, nClusts = j, 
        covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[3]])
    
    om_slog_hard_train <- cbind(om_slog_clustered_hard$train_data, as.data.frame(train_scores))
    om_slog_hard_test <- cbind(om_slog_clustered_hard$test_data, as.data.frame(test_scores))
      
    mse.seq_clusts[["om_slog_hard"]][i, j] <- fit_linear(
          train_data = om_slog_hard_train, 
          test_data = om_slog_hard_test)$mse

      
    # lcs hard
    lcs_clustered_hard <- hard_cluster(clusterward = clusterward_hard_lcs, nClusts = j, 
        covars=mvad_covars, y=num_month_em_last_year, train_idx, test_idx, dists[[2]])
    lcs_hard_train <- cbind(lcs_clustered_hard$train_data, as.data.frame(train_scores))
    lcs_hard_test <- cbind(lcs_clustered_hard$test_data, as.data.frame(test_scores))
      

    mse.seq_clusts[["lcs_hard"]][i, j] <- fit_linear(train_data = lcs_hard_train, test_data = lcs_hard_test)$mse
      
  # om-trate soft
    if (j < nSoftClusts) {
      om_trate_clustered_soft <- soft_cluster(dists[[1]], train_idx, test_idx = test_idx,
        nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
      
      if (om_trate_clustered_soft$converged) {
        om_trate_soft_train <- cbind(om_trate_clustered_soft$train_data, as.data.frame(train_scores))
        om_trate_soft_test <- cbind(om_trate_clustered_soft$test_data, as.data.frame(test_scores)) 
        mse.seq_clusts[["om_trate_soft"]][i,j] <- fit_linear(
            train_data = om_trate_soft_train, 
            test_data = om_trate_soft_test)$mse
        }
      
      # om-slog soft
      om_slog_clustered_soft <- soft_cluster(dists[[3]], train_idx, test_idx = test_idx,
        nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
      if (om_slog_clustered_soft$converged) {
        om_slog_soft_train <- cbind(om_slog_clustered_soft$train_data, as.data.frame(train_scores))
        om_slog_soft_test <- cbind(om_slog_clustered_soft$test_data, as.data.frame(test_scores))
        mse.seq_clusts[["om_slog_soft"]][i,j] <- fit_linear(
            train_data = om_slog_soft_train, 
            test_data = om_slog_soft_test)$mse
      }

      # lcs soft 
      lcs_clustered_soft <- soft_cluster(dists[[2]], train_idx, test_idx = test_idx,
        nClusts = j,fuzziness=fuzz_soft, covars = mvad_covars, y = num_month_em_last_year)
      if (lcs_clustered_soft$converged) {
        lcs_soft_train <- cbind(lcs_clustered_soft$train_data, as.data.frame(train_scores))
        lcs_soft_test <- cbind(lcs_clustered_soft$test_data, as.data.frame(test_scores))
        mse.seq_clusts[["lcs_soft"]][i,j] <- fit_linear(
            train_data = lcs_soft_train, 
            test_data = lcs_soft_test)$mse
      }
    }
  }
}


rmse.seqs_clusts <- list()
rmse.seqs_clusts$om_trate_hard <- sqrt(apply(mse.seq_clusts$om_trate_hard, 2, mean))
rmse.seqs_clusts$om_trate_soft <- sqrt(apply(mse.seq_clusts$om_trate_soft , 2, mean))
rmse.seqs_clusts$om_slog_soft <- sqrt(apply(mse.seq_clusts$om_slog_soft, 2, mean))
rmse.seqs_clusts$lcs_hard <- sqrt(apply(mse.seq_clusts$lcs_hard, 2, mean))
rmse.seqs_clusts$lcs_soft <- sqrt(apply(mse.seq_clusts$lcs_soft, 2, mean))



# Get best performance by each method (excpet sequence metrics because for that we 
# need best number of each type of clusters)
# square root after averaging
om_trate_hard_means_lm <- sqrt(apply(mse.cv.om_trate_hard, 2, mean))
om_trate_soft_means_lm <- sqrt(apply(mse.cv.om_trate_soft, 2, mean))
om_slog_hard_means_lm <- sqrt(apply(mse.cv.om_slog_hard, 2, mean))
om_slog_soft_means_lm <- sqrt(apply(mse.cv.om_slog_soft, 2, mean))
lcs_hard_means_lm <- sqrt(apply(mse.cv.lcs_hard, 2, mean))
lcs_soft_means_lm <- sqrt(apply(mse.cv.lcs_soft, 2, mean))
windows_means_lm <- sqrt(apply(mse.cv.windows, 2, mean))
harm_means_lm <- sqrt(apply(mse.cv.harm, 2, mean))

## Competition Evaluation
# Ultimately - we want the data and a plot 
rmse.frame <- data.frame(matrix(NA, nrow = 25, ncol = 8))
colnames(rmse.frame) <- c("OM T-Rate (Hard)", "OM T-Rate (Soft)", "OM INDELSLOG (Hard)", "OM INDELSLOG (Soft)", "LCS (Hard)", "LCS (Soft)", "Windows", "CFDA")

rmse.frame["OM T-Rate (Hard)"] <- om_trate_hard_means_lm
rmse.frame["OM T-Rate (Soft)"] <- om_trate_soft_means_lm
rmse.frame["OM INDELSLOG (Hard)"] <- om_slog_hard_means_lm
rmse.frame["OM INDELSLOG (Soft)"] <- om_slog_soft_means_lm
rmse.frame["LCS (Hard)"] <- lcs_hard_means_lm
rmse.frame["LCS (Soft)"] <- lcs_soft_means_lm
rmse.frame["Windows"] <- c(windows_means_lm, rep(NA,9))
rmse.frame["CFDA"] <- harm_means_lm
#rmse.frame["Sequence Metrics"] <- c(seq_mets_pc_means_lm, rep(NA, 13))
rmse.frame$Index <- 1:nrow(rmse.frame)

rmse_long <- rmse.frame |> filter(Index < 16) %>%
  pivot_longer(
    cols = -Index,          # Columns to pivot (all except Index)
    names_to = "Method",    # New column for the column names (OM, LCS, etc.)
    values_to = "RMSE"      # New column for the values
  )


pdf("LinearCompPlot.pdf",width=8,height=6)
ggplot(data = rmse_long, 
       aes(x = Index, y = RMSE, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  labs(title = "Linear Regression Performance by Method",
       x = "Components/Clusters",
       y = "RMSE Value (Averaged)",
       color = "Method") +
  theme_minimal()

dev.off()

saveRDS(mse.comp, file="LinearFullMSEs.rds")
saveRDS(cov.comp, file="LinearFullCoverages.rds")
saveRDS(mpiw.comp, file="LinearFullMPIW.rds")
saveRDS(mse.seq_clusts, file="LinearRMSEOptSeq+ClustMethods.rds")
saveRDS(rmse.frame, file="LinearCompetitionRMSE.rds")




  


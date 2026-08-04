# Functions for Sequence Outcome Prediction Work 

create_seq_data <- function(mvad) {
  mvad.labels <- c("employment", "further education", "higher education",
                  "joblessness", "school", "training")
  mvad.scodes <- c("EM","FE","HE","JL","SC","TR")

  # seqdef - creates a state sequence object (formulates sequences into labels)
  data.seq <- seqdef(mvad, 15:50, states=mvad.scodes, labels=mvad.labels)
  return(data.seq)

}

create_dists <- function(data.seq) {

  # creates distance matrices using three different sequence distance metrics. 
  # (1) OM-Trate (2) LCS (3) OM-INDELSLOG 
  
  # Calculate costs 
  costs.data <- seqcost(data.seq, method = "INDELSLOG", with.missing = TRUE)

  submat <- seqsubm(data.seq, method= "TRATE")
  # submat.2 <- submat
  # submat.2[] <- 1
  # diag(submat.2) <- 0
  
  # OM - T-RATE - Simply use the transition rates between each state. 
  dist.om_trate <- seqdist(data.seq, method="OM", indel = 1, sm = submat, with.missing = TRUE)
  
  # LCS
  dist.lcs <- seqdist(data.seq, method="LCS")
  
  # OM - SLOG (Must include costs)
  dist.om_slog <- seqdist(data.seq, method="OM", indel = costs.data$indel, sm = costs.data$sm, with.missing = TRUE)

  dists <- list(om_trate=dist.om_trate,lcs=dist.lcs, om_slog=dist.om_slog)
  
  return(dists)
}



create_rmetrics <- function(mvad.seq) {

  # encodings for the statebadness
  st_alphabet <- alphabet(mvad.seq)
  # Example: If alphabet is "A", "B", "C"
  st_prec_values <- c(1, -1, -1, 2, -1, -1) # A=1 (low badness), C=3 (high badness)
  names(st_prec_values) <- st_alphabet

  rmetrics <- data.frame(
      spells = seqindic(mvad.seq, "dlgth")$Dlgth,
      visited_states = seqindic(mvad.seq, "visited")$Visited,
      num_of_trans = seqindic(mvad.seq, "trans")$Trans,
      mean_spell_dur = seqindic(mvad.seq, "meand")$MeanD,
      # pedantic - this one pulled out a "seqivardur" "numeric" datatype (not sure
      # how it did both, so I have to force it to be numeric)
      sd_spell_dur = as.numeric(seqindic(mvad.seq, "dustd")$Dustd),
      # Diversity I
      entropy = seqindic(mvad.seq, "entr")$Entr,
      # more interested in states than spells (see paper)
      dss_subs = seqindic(mvad.seq, "nsubs")$Nsubs,
      complexity = seqindic(mvad.seq, "cplx")$Cplx,
      # could look at other turbulence measures
      turbulence = seqindic(mvad.seq, "turb")$Turb,
      badness = seqibad(seqdata = mvad.seq, stprec = st_prec_values),
      degradation = seqidegrad(seqdata = mvad.seq, stprec = st_prec_values),
      insecurity = seqinsecurity(seqdata = mvad.seq, stprec = st_prec_values)
    )
  
  return(rmetrics)
  
}


create_seqmets_pcs <- function(rmetrics, train_idx, test_idx, nSeqMetsPCs) {

  # Create Sequence Metric Principal Components
  train_rmetrics <- rmetrics[train_idx, ]
  test_rmetrics <- rmetrics[test_idx, ]

  pca_comps_train <- prcomp(x = train_rmetrics, center = TRUE, scale = TRUE)

  train_seq_scores <- pca_comps_train$x[, 1:nSeqMetsPCs]
  test_seq_scores <- predict(pca_comps_train, newdata = test_rmetrics)[,
  1:nSeqMetsPCs
  ]

  colnames(train_seq_scores) <- paste0("SeqMet", colnames(train_seq_scores))
  colnames(test_seq_scores) <- paste0("SeqMet", colnames(test_seq_scores))
  
  seq_data_sets <- list()
  seq_data_sets$train_pcs <- train_seq_scores
  seq_data_sets$test_pcs <- test_seq_scores

  return(seq_data_sets)

}


create_year_state_counts <- function(mvad) {

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


  return(mvad_states_wide)

}


create_bow_counts <- function(mvad) {

  # Make our windows of data - first three years
  mvad_states <- mvad[,c(1,15:50)]

  id_states <- mvad_states %>% pivot_longer(cols=colnames(mvad_states[2:37]), names_to="output") %>% count(id, value)


  mvad_states_wide <- id_states %>% 
    pivot_wider(id_cols=c(id), 
                names_from=value, values_from=n, values_fill = 0) %>% dplyr::select(-id)

  
  return(mvad_states_wide)

}

create_counts_pcs <- function(year_state_counts, mvad_covars, train_idx, test_idx, num_month_em_last_year) {


  pca_counts <- prcomp(x=year_state_counts[train_idx,], center=TRUE, scale=TRUE)
  train_scores <- pca_counts$x
  test_scores <- predict(pca_counts, newdata = year_state_counts[test_idx, ]) 

  colnames(train_scores) <- paste0("PC", 1:ncol(train_scores))
  colnames(test_scores) <- paste0("PC", 1:ncol(test_scores))

  counts_pcs <- list()
  counts_pcs$train_data <- cbind(mvad_covars[train_idx,] |> mutate(y = num_month_em_last_year[train_idx]), as.data.frame(train_scores))
  counts_pcs$test_data <- cbind(mvad_covars[test_idx,] |> mutate(y = num_month_em_last_year[test_idx]), as.data.frame(test_scores))

  return(counts_pcs)
}



create_long_data <- function(mvad) {

  M <- 36-1 #number of initial months-1

  mvad_wide <- mvad[,c(1,15:(15+M))] # Goes only to June 96 
  colnames(mvad_wide) <- c("id", as.character(0:M))

  # have time and id combinations (so basically only one state column)
  mvad_long <- gather(mvad_wide, key = time, value = state, "0":as.character(M))
  mvad_long$time <- as.numeric(mvad_long$time)
  mvad_long <- mvad_long %>% mutate(state=as.factor(state))

}


create_cfda_harms <- function(mvad_long, train_idx, test_idx, basis) {

  ids <- 1:712
  
  # compute encodings for - get warning messages that at least one states not in support of 
  # one basis function (I assume this is because of HE in the earlier years)
  fmca.train <- compute_optimal_encoding(mvad_long %>% filter(id %in% ids[train_idx]), basis, nCores = 7,verbose=F)
  pcs.train <- fmca.train$pc
  nComps <- ncol(pcs.train)
  colnames(pcs.train) <- paste0("PC",1:nComps)
  pcs.test <- predict(fmca.train,newdata=mvad_long %>% filter(id %in% ids[test_idx]),method="parallel",nCores=7)
  colnames(pcs.test) <- paste0("PC",1:nComps)
  
  harm_data <- list()
  harm_data$train_data <- mvad_covars[train_idx,]%>% mutate(y=num_month_em_last_year[train_idx]) %>% add_column(as_tibble(pcs.train)) 
  harm_data$test_data <- mvad_covars[test_idx,] %>% mutate(y=num_month_em_last_year[test_idx]) %>% add_column(as_tibble(pcs.test))

  return(harm_data)
}


# Assigns test data to a hard clustering solution based on a passed in distance matrix
assign_new <- function(dist_mat,train_idx,train.clust,test_idx) {
  
  # Empty cluster assignments
  len <- length(train_idx) 
  clusts <- rep(0,len) 
  
  # Add training assignments assigned 
  clusts[train_idx] <- train.clust 
  
  # calculate distance to each cluster 
  dist2clust <- apply(dist_mat,1,function(x,cl) {tapply(x,cl,mean)},cl=clusts)  
  newClusts <- apply(dist2clust[-1,],2,which.min) # drop first row (unassigned cluster (0))
  
  return(newClusts[test_idx])
  
}

# Creates test set membership probabilities based on k nearest neighbors in the training set. 
assign_new_fanny_knn <- function(dist_mat,train.idx,train.memb,test.idx,memb.exp=1.75, k = 10) {

  len <- length(train.idx) # all members
  test_locs <- (1:len)[test.idx] # indices of test locations

  K <- ncol(train.memb) # number of clusters
  M <- length(test_locs) # number of new obs to assign
  test.memb <- matrix(NA,M,K) # empty membership probability matrix

  # Exponent for distances (Weight of membership expression)
  dist_exp <- -2 / (memb.exp - 1)

   for (i in seq_along(test_locs)) {
    
    test_loc <- test_locs[i] # get current test location in mvad 
    dist_to_test <- dist_mat[test_loc, train.idx] # get distances from training points to test points 
    train_nn <- dist_to_test[order(dist_to_test)[1:k]] # indexes in the training set (not mvad)
     
    train_nn_idx <- names(train_nn) # use names because the train.memb is indexed by mvad (712)
    train_nn_dists <- as.vector(train_nn)
     
    if (any(train_nn_dists == 0)) {

      # If any neighbors have a distance of 0, assign that point's membership probailities to the test point
      nearestn <- train_nn_idx[1] # get the name of the original index in mvad 
      test.memb[i,] <- train.memb[nearestn, ]
    } else {

      test_dists_term <- train_nn_dists ^ dist_exp
      nn_memb_probs <- train.memb[train_nn_idx,,drop=FALSE]
      nn_membs <- colSums(sweep(nn_memb_probs, MARGIN = 1, STATS = test_dists_term, FUN = "*")) / sum(test_dists_term)
      
      test.memb[i,] <- nn_membs
    }
   }

  return(test.memb)
}

# Function to fit random forest and calculate statistics 
fit_rf <- function(train_data, test_data, mtry = NA, num.trees = 1000) {

  if (!is.na(mtry)) {
    mtry_tr = mtry
  } else {
    mtry_tr = ncol(train_data) - 1
  }
 
  # fit with hard clusters 
  fit.rf <- ranger(y ~ .,
                        data = train_data,
                        num.trees = num.trees, 
                        mtry = mtry_tr,
                        keep.inbag = TRUE,
                        respect.unordered.factors = "order",
                        classification = FALSE,
                        quantreg = TRUE)


  pred.rf <- predict(fit.rf, data = test_data)
  pred_quantiles.rf <- predict(fit.rf, data = test_data, type = "quantiles", quantiles = c(0.025, 0.975))
  mse.rf <- calc_mse(preds = pred.rf$predictions, test_data$y)

  # pull out quantiles 
  lwr <- pred_quantiles.rf$predictions[, 1]
  upr <- pred_quantiles.rf$predictions[, 2]

  mpiw <- mean(upr - lwr)
  coverage <- mean(test_data$y >= lwr & test_data$y <= upr)   

  return(list(mse=mse.rf, mpiw=mpiw, coverage=coverage))
  
}


# Perform a hard clustering  using the training data, and assign test data to clusters based on distances
# Returns data with no covariates 
hard_cluster <- function(clusterward, nClusts, covars = NULL, y, train_idx, test_idx, dist_matrix) {
  # Takes in a clusterward (agnes) clustering and cuts it into a number of clusters (nClust)
  # dist_matrix is contains the distances between all data points


  cut1 <- cutree(clusterward,k=nClusts)
  clust1.fac <- factor(cut1)
  train_data <-  data.frame(cluster=clust1.fac, y=y[train_idx])
  
  cut2 <- assign_new(dist_matrix,train_idx,cut1,test_idx)
  clust2.fac <- factor(cut2)
  test_data <- data.frame(cluster=clust2.fac, y=y[test_idx])

  if (!is.null(covars)) {
    train_data <- cbind(covars[train_idx, ], train_data)
    test_data <- cbind(covars[test_idx,], test_data)

  }
  
  return(list(train_data=train_data, test_data=test_data))
}


# soft_cluster_safeno <- function(dist_matrix,train_idx,test_idx, nClusts, fuzziness=1.75, covars, y, prior.n.equiv=100, max_iters = 500, tol=1e-15) {
  
#   clustering_soft <- fanny(dist_matrix[train_idx, train_idx], 
#                            k=nClusts, memb.exp=fuzziness, diss=TRUE, maxit = max_iters, tol=tol)
  
#   train.memb <- clustering_soft$membership
#   test.memb <- assign_new_fanny(dist_mat=dist_matrix,train.idx=train_idx, train.memb=train.memb, 
#                                test.idx=test_idx,memb.exp=fuzziness,n.draws=200,
#                                prior.n.equiv=prior.n.equiv) 

#   converged <- (clustering_soft$convergence["converged"] == 1)

#   safe_condition_no <- condition_no_check(clustering_soft$membership)$safe_condition_no
#     condition_no <- condition_no_check(clustering_soft$membership)$condition_no
  
#   colnames(train.memb) <- paste0("Cluster",1:nClusts)
#   colnames(test.memb) <- paste0("Cluster",1:nClusts)
  
#   if (nClusts == 2) {
#     train_data <- covars[train_idx,] %>% 
#       mutate(y=y[train_idx], cluster1=train.memb[,-1])
#     test_data <- covars[test_idx,] %>% 
#       mutate(y=y[test_idx], cluster1=test.memb[,-1])
#   } else {
#     train_data <- cbind(covars[train_idx,], train.memb[,-1]) %>%
#       mutate(y=y[train_idx])
#     test_data <- cbind(covars[test_idx,], test.memb[,-1]) %>%
#       mutate(y=y[test_idx])
#   } 
#   return(list(train_data=train_data, test_data=test_data, converged=converged, safe_condition_no=safe_condition_no, condition_no=condition_no))
# }


soft_cluster <- function(dist_matrix,train_idx,test_idx, nClusts, fuzziness=1.75, covars, y, max_iters = 2000, tol=1e-10, knn_soft_assign=15) {
  

  clustering_soft <- fanny(dist_matrix[train_idx, train_idx], 
                           k=nClusts, memb.exp=fuzziness, diss=TRUE, maxit = max_iters, tol=tol)
  
  train.memb <- clustering_soft$membership
  test.memb <- assign_new_fanny_knn(dist_mat=dist_matrix,train.idx=train_idx, train.memb=train.memb, 
                               test.idx=test_idx,memb.exp=fuzziness,
                               k=knn_soft_assign) 



  converged <- (clustering_soft$convergence["converged"] == 1)
  
  colnames(train.memb) <- paste0("Cluster",1:nClusts)
  colnames(test.memb) <- paste0("Cluster",1:nClusts)
  
  if (nClusts == 2) {
    train_data <- covars[train_idx,] %>% 
      mutate(y=y[train_idx], cluster1=train.memb[,-1])
    test_data <- covars[test_idx,] %>% 
      mutate(y=y[test_idx], cluster1=test.memb[,-1])
  } else {
    train_data <- cbind(covars[train_idx,], train.memb[,-1]) %>%
      mutate(y=y[train_idx])
    test_data <- cbind(covars[test_idx,], test.memb[,-1]) %>%
      mutate(y=y[test_idx])
  } 
  return(list(train_data=train_data, test_data=test_data, converged=converged)) 
}


# Check the condition number of assigned membership probabilities matrix 
condition_no_check <- function(memb_probs, threshold = 100) {

  safe_condition_no <- TRUE

  condition_no <- kappa(memb_probs, exact=TRUE)

  if (condition_no > threshold) {
    safe_condition_no <- FALSE
  }

  condition_no_info <- list(condition_no = condition_no, safe_condition_no = safe_condition_no)

  return(condition_no_info)
}


soft_cluster_sim <- function(dist_matrix,train_idx,test_idx, nClusts, fuzziness=1.75, y, max_iters = 2000, tol=1e-10, knn_soft_assign=15) {
  
  clustering_soft <- fanny(dist_matrix[train_idx, train_idx], 
                           k=nClusts, memb.exp=fuzziness, diss=TRUE, maxit = max_iters, tol=tol)
  
  train.memb <- clustering_soft$membership
  test.memb <- assign_new_fanny_knn(dist_mat=dist_matrix,train.idx=train_idx, train.memb=train.memb, 
                               test.idx=test_idx,memb.exp=fuzziness,
                               k=knn_soft_assign) 

  converged <- (clustering_soft$convergence["converged"] == 1)
  condition_no_output <- condition_no_check(clustering_soft$membership)

  safe_condition_no <- condition_no_output$safe_condition_no
  condition_no <- condition_no_output$condition_no
  
  colnames(train.memb) <- paste0("Cluster",1:nClusts)
  colnames(test.memb) <- paste0("Cluster",1:nClusts)
  
  if (nClusts == 2) {
    train_data <- data.frame(Cluster1=train.memb[,-1], y=y[train_idx])
    test_data <- data.frame(Cluster1=test.memb[,-1], y=y[test_idx])
  } else {
    train_data <- data.frame(train.memb[,-1]) |> mutate(y= y[train_idx])
    test_data <- data.frame(test.memb[,-1]) |> mutate(y= y[test_idx])
  } 

  return(list(train_data=train_data, test_data=test_data, converged=converged, safe_condition_no=safe_condition_no, condition_no = condition_no))
}


# Calculate the mean squared error 
calc_mse <- function(preds, y_dat) {
  return(sum((preds - y_dat)^2)/length(y_dat)) 
}

# Fits a linear model and calculates mse, mean interval prediction width and coverage
fit_linear <- function(train_data, test_data) {

  fit.lin <- lm(y~., data=train_data) 
  preds.lin <- predict(fit.lin, newdata = test_data, interval = "prediction", level = 0.95)
  mse.lin <- calc_mse(preds.lin[,"fit"], y_dat=test_data[,"y"])
  
  # Intervals 
  lwr <- preds.lin[, "lwr"]
  upr <- preds.lin[, "upr"]
  mpiw <- mean(upr - lwr)

  coverage <- mean((test_data$y >= lwr) & (test_data$y <= upr))

  return(list(mse=mse.lin, mpiw=mpiw, coverage=coverage))
  
}

# Function for calculating minimum (with NA values)
safe_min <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(min(x, na.rm = na.rm))
  }
}

# Function for calculating which minimum (with NA values)
safe_which_min <- function(x) {
  if (all(is.na(x))) {
    return(NA) # Returns NA instead of integer(0)
  } else {
    return(which.min(x))
  }
}

# Function for calculating mode (with missing values)
safe_mode <- function(x) {
  if (all(is.na(x))) {
    return(NA) # Returns NA instead of integer(0)
  } else {
    if (length(unique(x)) == length(x)) {
      return(round(mean(x), 0))
    } else {
      tab <- table(x)
      return(as.integer(names(tab)[which.max(tab)]))
    }
  }

}





# Plots helper functions

# MSE Plot
get_method_metrics <- function(data_list, method_name, method_labs) {

  # pull out cross validated mses for given method
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)

  method_cov <- lapply(data_list$cov, function(x) x[[method_name]])
  combined_cov <- simplify2array(method_cov)
  
  method_mpiw <- lapply(data_list$mpiw, function(x) x[[method_name]])
  combined_mpiw <- simplify2array(method_mpiw)



  if (method_name == "demos") {
    mean_mse <- mean(combined_mse)
    mean_cov <- mean(combined_cov)
    mean_mpiw <- mean(combined_mpiw)

  } else {
    # find means across the folds (of all cross )
    mean_mse <- apply(combined_mse, 2, mean, na.rm = TRUE)
    mean_cov <- apply(combined_cov, 2, mean, na.rm = TRUE)
    mean_mpiw <- apply(combined_mpiw, 2, mean, na.rm = TRUE)
  }


  if (method_name %in% c(
  "om_trate_soft", 
  "om_slog_soft",
  "lcs_soft"
)) {
  method_con <- lapply(data_list$convergence, function(x) x[[method_name]])
  combined_con <- simplify2array(method_con)
  
  # Number of Convergences 
  convergences <- apply(combined_con, 2, sum)
  
  } else {
    convergences <- rep(NA, length(mean_mse))
  }

  

  df <- data.frame(
      Method = method_name, 
      Index = 1:length(mean_mse),
      RMSE   = sqrt(mean_mse), 
      COV    = mean_cov, 
      MPIW   = mean_mpiw,
      Convergences = convergences, 
      stringsAsFactors = FALSE)


  df <- df |> mutate(Method = recode(Method, !!!method_labs))

  return(df)
}




# function to get best mse (considering mtry) for every method/component combination
get_nonlinear_metrics_methods <- function(data_list, method_name) {
  
  # Extract the list of 20 CV 3d arrays (fold, comps, mtry)
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)
  method_cov <- lapply(data_list$cov, function(x) x[[method_name]])
  combined_cov <- simplify2array(method_cov)
  method_mpiw <- lapply(data_list$mpiw, function(x) x[[method_name]])
  combined_mpiw <- simplify2array(method_mpiw)
  
  # within every cross validation - find mse average by mtry 
  avg_cv_mse <- apply(combined_mse, c(2, 3, 4), mean, na.rm = TRUE)
  avg_cv_cov <- apply(combined_cov, c(2, 3, 4), mean, na.rm = TRUE)
  avg_cv_mpiw <- apply(combined_mpiw, c(2, 3, 4), mean, na.rm = TRUE)

  # essentially the location maps to the mtry values
  best_mtry_cv <- apply(avg_cv_mse, c(1, 3), safe_which_min)
  modal_mtrys <- apply(best_mtry_cv, 1, safe_mode)

  # the actual min corresponds to this so we will get mses that correspond to the best
  best_avg_mse_comps_cv <- apply(avg_cv_mse, c(1, 3), safe_min)
  # now we can average these mses (picked by mtry) to find mse by comp
  rmse_comps <- sqrt(apply(best_avg_mse_comps_cv, 1, mean, na.rm=TRUE))

  # data frame that has the indexes needed to get the mtrys for each components/cv combination
  best_mtry_indices <- cbind(
    as.vector(row(best_mtry_cv)), 
    as.vector(best_mtry_cv), 
    as.vector(col(best_mtry_cv))
  )

  # copy dimensions of avg_cv_cov 
  dims <- dim(avg_cv_cov)

if (method_name %in% c("om_trate_soft", "om_slog_soft", "lcs_soft")) {
    method_con <- lapply(data_list$convergence, function(x) x[[method_name]])
    combined_con <- simplify2array(method_con)
    
    convergences <- apply(combined_con, 2, sum) 
  
} else {
    convergences <- rep(NA, dims[1])
  }

  best_avg_cov_comps_cv <- matrix(avg_cv_cov[best_mtry_indices], nrow = dims[1], ncol = dims[3])
  best_avg_mpiw_comps_cv <- matrix(avg_cv_mpiw[best_mtry_indices], nrow = dims[1], ncol = dims[3])

  cov_comps <- apply(best_avg_cov_comps_cv, 1, mean)
  mpiw_comps <- apply(best_avg_mpiw_comps_cv, 1, mean)

  return(list(rmses = rmse_comps, cov = cov_comps, mpiw =  mpiw_comps, modal_mtry = modal_mtrys, mtry_dis = best_mtry_cv, 
    method_name = method_name, convergences=convergences)) 

}


# Separate function for demographics (as the number of components does not change)
get_nonlinear_demos <- function(data_list, method_name) {
  
  # Extract the list of 20 CV 3d arrays (fold, comps, mtry)
  method_mse <- lapply(data_list$mse, function(x) x[[method_name]])
  combined_mse <- simplify2array(method_mse)
  method_cov <- lapply(data_list$cov, function(x) x[[method_name]])
  combined_cov <- simplify2array(method_cov)
  method_mpiw <- lapply(data_list$mpiw, function(x) x[[method_name]])
  combined_mpiw <- simplify2array(method_mpiw)
  
  # within every cross validation - find mse average by mtry 
  avg_cv_mse <- apply(combined_mse, c(2, 3), mean, na.rm = TRUE)
  avg_cv_cov <- apply(combined_cov, c(2, 3), mean, na.rm = TRUE)
  avg_cv_mpiw <- apply(combined_mpiw, c(2, 3), mean, na.rm = TRUE)

  # essentially the location maps to the mtry values
  best_mtry_cv <- apply(avg_cv_mse, 2, safe_which_min)
  modal_mtry <- safe_mode(best_mtry_cv)

  # the actual min corresponds to this so we will get mses that correspond to the best
  best_avg_mse_comps_cv <- apply(avg_cv_mse, 2, safe_min)
  # now we can average these mses (picked by mtry) to find mse by comp
  rmse_comps <- sqrt(mean(best_avg_mse_comps_cv))

  cov_comps <- mean(avg_cv_cov[best_mtry_cv])
  mpiw_comps <- mean(avg_cv_mpiw[best_mtry_cv])

  df <- data.frame(comps=11, rmse=rmse_comps, cov=cov_comps, mpiw=mpiw_comps, method="Demographics", modal_mtry=modal_mtry, convergences=NA)
  df$mtry_dist <- list(best_mtry_cv)
  return(df)

}


get_nonlinear_metrics <- function(data_list, method_labels) {

  # Create dataframe of nonlinear competition 
  methods_metrics <- map_dfr(names(method_labels), function(m) {
    if (m != "demos") {

    res <- get_nonlinear_metrics_methods(data_list, m)  
    data.frame(
      comps = seq_along(res$rmses),
      rmse = res$rmses,
      cov = res$cov, 
      mpiw = res$mpiw,
      method = method_labels[[m]], 
      modal_mtry = res$modal_mtry,
      mtry_dist = I(split(res$mtry_dis, row(res$mtry_dis))),
      convergences = res$convergences
    )
  }
    })
  
  
  demos_metrics <- get_nonlinear_demos(data_list, "demos")

  nonlin_metric_data <- rbind(methods_metrics, demos_metrics)
  rownames(nonlin_metric_data) <- 1:nrow(nonlin_metric_data)

  return(nonlin_metric_data)
  
}



var_imp_plot <- function(mvad_data, best_mtry, title_string, subtitle_string = "", custom_labs = FALSE) {

  fit.rf <- ranger(y ~ .,
                        data = mvad_data,
                        num.trees = 1000, 
                        mtry = best_mtry,
                        keep.inbag = TRUE,
                        respect.unordered.factors = "order",
                        classification = FALSE,
                        importance = "impurity",
                        quantreg = TRUE)
  
  var_imps <- data.frame(fit.rf$variable.importance) |> rownames_to_column("covar") |> 
    rename(var_imp = fit.rf.variable.importance)

  p <- ggplot(data = var_imps, aes(x=var_imp, y=reorder(covar, var_imp))) + geom_col(fill="steelblue") + 
    labs(x="Variable Importance", y="Covariate", title=title_string, subtitle=subtitle_string)

  return(p)
}




vivi_variable_calcs <- function(mvad_data, best_mtry) {
  
    fit.rf <- ranger(y ~ .,
                        data = mvad_data,
                        num.trees = 1000, 
                        mtry = best_mtry,
                        keep.inbag = TRUE,
                        respect.unordered.factors = "order",
                        classification = FALSE,
                        importance = "impurity",
                        quantreg = TRUE)
  
  
  viviRangEmbedded <- vivi(fit = fit.rf , 
                         data = mvad_data, 
                         response = "y",
                         importanceType = "impurity")
  
  
  return(viviRangEmbedded)
  
}
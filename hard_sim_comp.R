# Simulated Data Evaluation of Life Course Methods

library(tidyverse)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(cfda)
library(foreach)
library(doParallel)

source("mvad_seqout_functions.R")


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
nComps <- 8
nSoft <- 5
nSets <- 150
fuzz_soft <- 1.25

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

# locations for training and testing sets - first 900 
# are training, second 900 are testing. 
set_train_idx <- c(rep(TRUE, 900), rep(FALSE, 900))
set_test_idx <- !set_train_idx

# initialize basis functions - can do this generally
basis <- create.bspline.basis(c(1, 12), nbasis = 3, norder = 3)

# data structure for holding the best number of 
# components by performance method. 
best_n_comps <- list()
best_n_comps$concord <- array(NA, c(nMethods, nSets))
best_n_comps$indep_semi <- array(NA, c(nMethods, nSets))
best_n_comps$indep <- array(NA, c(nMethods, nSets))
best_n_comps$semi_concord <- array(NA, c(nMethods, nSets))
best_n_comps$semi <- array(NA, c(nMethods, nSets))


set_conv_train_track <- list()
set_conv_train_track$om_trate <- set_conv_train_track$om_slog <- set_conv_train_track$lcs <- array(FALSE, c(folds, nSoft))

train_conv_tracker <- list()
train_conv_tracker$concord <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$indep_semi <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$indep <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$semi_concord <- replicate(nSets, set_conv_train_track, simplify = FALSE)
train_conv_tracker$semi <- replicate(nSets, set_conv_train_track, simplify = FALSE)


for (file_name in file_names) {
  print(file_name)
  # Load in the file 
  file_path <- paste0("sim_data/", sim_type, "/", file_name)
  load(file_path)
  cat("Loaded:", file_path, "\n")

  sim_group <- sub("^[^_]*_(.*)\\..*$", "\\1", file_name)
  
  results <- foreach(n = 1:nSets, .packages = c("tidyverse", "cluster", "TraMineR", "cfda", "TraMineRextras")) %dopar% {
    
    cat("\nProcessing", file_path, "dataset", n)
    
    # pick set to do cross validation on 
    data_wide <- df_long_list[[n]] %>% group_by(ID) %>% 
      mutate(y = y[time=="T1"]) %>% ungroup() %>%  # make every y value equal to the "first" value  
      pivot_wider(id_cols=c(ID, Class, y), names_from=time, values_from=state) %>% 
      dplyr::select(-Class, -ID) # drop Class variable and ID is the same as the index 
    
    # create distance matrices (first column is Y value) 
    dists <- create_dists(data.seq=seqdef(data_wide, 2:13, 
                      states=c("A", "B", "C")))
    
    set_train_data <- data_wide[set_train_idx,]
    set_test_data <- data_wide[set_test_idx,]
    
  
    # Re-initialize for every cv (every set)
    mse.cv <- array(NA, c(folds, nComps, nMethods))
    
    # CFDA prep scetion ƒ
    cfda_data_long <- df_long_list[[n]] %>% mutate(
      time= as.numeric(str_sub(time, start = 2, end = -1)),
      state=as.factor(state)) %>%
      rename(id=ID) %>%
      dplyr::select(-Class, -y)

    cfda_long_set_train <- cfda_data_long %>% filter(id %in% 1:900)
    cfda_long_set_test <- cfda_data_long %>% filter(id %in% 901:1800)


    # count across dataset
    set_data_windows <- cfda_data_long %>% mutate(
      window = if_else(time < 7, 1, 2)) %>%
      count(id, window, state) %>%
      pivot_wider(
        id_cols=c(id), names_from=c(state,window), values_from=n, values_fill = 0)

    
    # soft clustering convergence tracking 
    set_conv_train_track <- list()
    set_conv_train_track$om_trate <- set_conv_train_track$om_slog <- set_conv_train_track$lcs <- array(FALSE, c(folds, nSoft))

    
    # cross validation
    for(i in 1:folds) {
        cat("Fold Number ",i,"\n")
        cv_test_idx <- cv_idx == i
        cv_train_idx <- !cv_test_idx
        
        # get fold train and test sets, from training set.
        y.cv_train <- set_train_data[cv_train_idx, "y"]
        y.cv_test <- set_train_data[cv_test_idx, "y"]
        
        # hard clustering - will get "cut" depending on number of clusters. 
        clusterward_trate_hard <- agnes(
          dists[[1]][set_train_idx, set_train_idx][cv_train_idx,cv_train_idx],
            diss=TRUE, method="ward")

              # hard clustering - will get "cut" depending on number of clusters. 
        clusterward_slog_hard <- agnes(
          dists[[3]][set_train_idx, set_train_idx][cv_train_idx,cv_train_idx],
            diss=TRUE, method="ward")
      
              # hard clustering - will get "cut" depending on number of clusters. 
        clusterward_lcs_hard <- agnes(
          dists[[2]][set_train_idx, set_train_idx][cv_train_idx,cv_train_idx],
            diss=TRUE, method="ward")

  
        # one basis function (I assume this is because of HE in the earlier years)
        fmca.train <- compute_optimal_encoding(cfda_long_set_train %>% filter(id %in% seq(1:900)[cv_train_idx]),
                                              basis,verbose=F)
        pcs.train <- fmca.train$pc
        nCFDAComps <- ncol(pcs.train)

        pcs.test <- predict(fmca.train,newdata=
                              cfda_long_set_train %>% filter(id %in% seq(1:900)[cv_test_idx]),verbose=F)

        #it's relatively easy to run this with 5 or 10 harmonics - or in the equation, below
        cv_train_harm <- data.frame(y = y.cv_train) %>% add_column(as_tibble(pcs.train))
        cv_test_harm <- data.frame(y = y.cv_test) %>% add_column(as_tibble(pcs.test))
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

        # 5 folds and 5 potential number of solutions                      
        soft_convergenece_train_tracker <- array(NA, c(folds, 5))
      
        for(j in 1:nComps) {
          
          if (j > 1) {
            # ***************
            # hard clustering. 
            om_trate_clusters <- hard_cluster_sim(clusterward_trate_hard, j, set_train_data[["y"]], 
              cv_train_idx, cv_test_idx, dist_matrix = dists[[1]][set_train_idx, set_train_idx])

            om_slog_clusters <- hard_cluster_sim(clusterward_slog_hard, j, set_train_data[["y"]], 
              cv_train_idx, cv_test_idx, dist_matrix = dists[[3]][set_train_idx, set_train_idx])
            
            lcs_clusters <- hard_cluster_sim(clusterward_lcs_hard, j, set_train_data[["y"]], 
              cv_train_idx, cv_test_idx, dist_matrix = dists[[2]][set_train_idx, set_train_idx])


            # hard is method "1" - fit model and calculate mse
            mse.cv[i, j, 3] <- fit_linear(
              train_data = om_trate_clusters$train_data,
              test_data = om_trate_clusters$test_data)
            mse.cv[i, j, 4] <- fit_linear(
              train_data = om_slog_clusters$train_data,
              test_data = om_slog_clusters$test_data)
            mse.cv[i, j, 5] <- fit_linear(
              train_data = lcs_clusters$train_data,
              test_data = lcs_clusters$test_data)
            
          # **********************
          # soft clustering - we have to consider potential convergence problems
            if (j <= nSoft) {
              om_trate_soft_train <- soft_cluster_sim(dists[[1]][set_train_idx, set_train_idx],
                cv_train_idx, cv_test_idx, nClusts=j, fuzziness = fuzz_soft, y=set_train_data[["y"]]
              )
              if (om_trate_soft_train$converged) {
                # hard is method "2" - fit model and calculate mse 
                mse.cv[i, j, 6] <- fit_linear(om_trate_soft_train$train_data, om_trate_soft_train$test_data)
                set_conv_train_track[["om_trate"]][i,j] <- TRUE
              } else {
                mse.cv[i, j, 6] <- NA
              }

              om_slog_soft_train <- soft_cluster_sim(dists[[3]][set_train_idx, set_train_idx],
                cv_train_idx, cv_test_idx, nClusts=j, fuzziness = fuzz_soft, y=set_train_data[["y"]]
              )

              if (om_slog_soft_train$converged) {
                # hard is method "2" - fit model and calculate mse 
                mse.cv[i, j, 7] <- fit_linear(om_slog_soft_train$train_data, om_slog_soft_train$test_data)
                set_conv_train_track[["om_slog"]][i,j] <- TRUE

              } else {
                mse.cv[i, j, 7] <- NA
              }

              lcs_soft_train <- soft_cluster_sim(dists[[2]][set_train_idx, set_train_idx],
                cv_train_idx, cv_test_idx, nClusts=j, fuzziness = fuzz_soft, y=set_train_data[["y"]]
              )

              if (lcs_soft_train$converged) {
                # hard is method "2" - fit model and calculate mse 
                mse.cv[i, j, 8] <- fit_linear(lcs_soft_train$train_data, lcs_soft_train$test_data)
                set_conv_train_track[["lcs"]][i,j] <- TRUE

              } else {
                mse.cv[i, j, 8] <- NA
              }
            }
          } 
        
          # ****************
          # CFDA 
          
          # y value is the first column, so we index to j + 1 to get j elements for this round.
          # j = 2 means we look at the first PC (so clusters we get 2 but PCs we get 1 for j=2)
          # "3" is CFDA
          mse.cv[i, j, 1] <- fit_linear(cv_train_harm[,1:(j+1)], cv_test_harm[,1:(j+1)])

          # # ****************
          # # Windows

          # # can't do more than 6 PCs because there are only 6 dimensions
          # # (should really do less than 6)
          if (j < 7) {
            mse.cv[i, j, 2] <- fit_linear(cv_train_windows[,1:(j+1)], cv_test_windows[,1:(j+1)])
          }
        }
    }
    
    print("folds done")
    
    # Fitting by best parameters
    best_cfda_harms <- which.min(apply(mse.cv[,,1], 2, mean, na.rm = TRUE))
    best_windows_pcs <- which.min(apply(mse.cv[,,2], 2, mean, na.rm = TRUE))
    best_trate_hard_clusts <- which.min(apply(mse.cv[,,3], 2, mean, na.rm = TRUE))
    best_slog_hard_clusts <- which.min(apply(mse.cv[,,4], 2, mean, na.rm = TRUE))
    best_lcs_hard_clusts <- which.min(apply(mse.cv[,,5], 2, mean, na.rm = TRUE))
    
    # need a little more nuance here in terms of clustering 
    # if there are any NAs in a column (ie in any of the folds) we do not look 
    # at that number of clusters (i.e. we ignore it and look at the average for 
    # number of components that have converged). 
    best_trate_soft_clusts <- which.min(apply(
      data.frame(mse.cv[,,6]) |> mutate(across(where(~any(is.na(.x))), ~NA)), 2, mean, na.rm = TRUE))[[1]]
    best_slog_soft_clusts <- which.min(apply(
      data.frame(mse.cv[,,7]) |> mutate(across(where(~any(is.na(.x))), ~NA)), 2, mean, na.rm = TRUE))[[1]]
    best_lcs_soft_clusts <- which.min(apply(
      data.frame(mse.cv[,,8]) |> mutate(across(where(~any(is.na(.x))), ~NA)), 2, mean, na.rm = TRUE))[[1]]

    
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
      dists[[3]][set_train_idx, set_train_idx],
      diss=TRUE, method="ward")
  
    om_trate_clusters_set <- hard_cluster_sim(clusterward_om_trate_hard_set,
      best_trate_hard_clusts, data_wide[["y"]], set_train_idx, set_test_idx, dists[[1]])
    om_slog_clusters_set <- hard_cluster_sim(clusterward_om_slog_hard_set,
      best_slog_hard_clusts, data_wide[["y"]], set_train_idx, set_test_idx, dists[[3]])
    lcs_clusters_set <- hard_cluster_sim(clusterward_lcs_hard_set,
      best_lcs_hard_clusts, data_wide[["y"]], set_train_idx, set_test_idx, dists[[2]])
  
    mse_om_trate_hard <- fit_linear(
      train_data=om_trate_clusters_set$train_data, 
      test_data=om_trate_clusters_set$test_data)
    mse_om_slog_hard <- fit_linear(
      train_data=om_slog_clusters_set$train_data, 
      test_data=om_slog_clusters_set$test_data)
    mse_lcs_hard <- fit_linear(
      train_data=lcs_clusters_set$train_data, 
      test_data=lcs_clusters_set$test_data)
    
    
    # om trate soft 
    # if there is no training convergence, mse is NA 
    if (!identical(best_trate_soft_clusts, integer(0))) {
      om_trate_soft_set <- soft_cluster_sim(dists[[1]],
              set_train_idx, set_test_idx, nClusts=best_trate_soft_clusts, fuzziness = fuzz_soft, y=data_wide[["y"]])
      if (om_trate_soft_set$converged) {
        mse_om_trate_soft <- fit_linear(train_data=om_trate_soft_set$train, test_data=om_trate_soft_set$test_data)
      }
    } else {
      mse_om_trate_soft <- NA
    }

    # om-slog soft 
    if (!identical(best_slog_soft_clusts, integer(0))) {
      om_slog_soft_set <- soft_cluster_sim(dists[[3]],
              set_train_idx, set_test_idx, nClusts=best_slog_soft_clusts, fuzziness = fuzz_soft, y=data_wide[["y"]])
      if (om_slog_soft_set$converged) {
        mse_om_slog_soft <- fit_linear(train_data=om_slog_soft_set$train, test_data=om_slog_soft_set$test_data)
      }
    } else {
      mse_om_slog_soft <- NA
    }


    # lcs soft
    if (!identical(best_trate_soft_clusts, integer(0))) {
      lcs_soft_set <- soft_cluster_sim(dists[[2]],
              set_train_idx, set_test_idx, nClusts=best_lcs_soft_clusts, fuzziness = fuzz_soft, y=data_wide[["y"]])
      if (lcs_soft_set$converged) {
        mse_lcs_soft <- fit_linear(train_data=lcs_soft_set$train, test_data=lcs_soft_set$test_data)
      }
    } else {
      mse_lcs_soft <- NA
    }


    # CFDA - refit
    fmca.train_set <- compute_optimal_encoding(cfda_long_set_train, basis,verbose=F)
    pcs.train_set <- fmca.train_set$pc
    pcs.test_set <- predict(fmca.train_set,newdata=cfda_long_set_test,verbose=F)

    set_train_harm <- data.frame(y = y.train_set) %>% add_column(as_tibble(pcs.train_set[,1:best_cfda_harms]))
    set_test_harm <- data.frame(y = y.test_set) %>% add_column(as_tibble(pcs.test_set[,1:best_cfda_harms]))
    colnames(set_test_harm) <- colnames(set_train_harm) # ensure names match

    mse_cfda <- fit_linear(set_train_harm, set_test_harm)
    

    # Windows - refit
    pca_windows_train_set <- prcomp(x=set_data_windows %>% filter(id %in% 1:900) %>% dplyr::select(-id))
    set_train_scores <- pca_windows_train_set$x
    set_test_scores <- predict(pca_windows_train_set,
                              newdata = set_data_windows %>% filter(id %in% 901:1800) %>% dplyr::select(-id))
    
    # add principal components to demographic data (no longer need original windows)
    set_train_windows <- cbind(as.data.frame(y.train_set), set_train_scores)
    set_test_windows <- cbind(as.data.frame(y.test_set), set_test_scores)
    
    mse_wind <- fit_linear(
      train_data=cbind(as.data.frame(y.train_set), set_train_scores),
      test_data=cbind(as.data.frame(y.test_set), set_test_scores))   
    
  list(
        mses = c(mse_cfda, mse_wind, mse_om_trate_hard, mse_om_slog_hard, mse_lcs_hard,mse_om_trate_soft, 
          mse_om_slog_soft, mse_lcs_soft),
        best_k = c(best_cfda_harms, best_windows_pcs, best_trate_hard_clusts, best_slog_hard_clusts,
          best_lcs_hard_clusts, best_trate_soft_clusts,best_slog_soft_clusts, best_lcs_soft_clusts), 
        set_conv_train_track = set_conv_train_track
      )
  }
  # Transfer results from workers back to global matrices
  # 5. Populate global matrices from results list

  for(n in 1:nSets) {
    mse.sims[[sim_group]][,n] <- results[[n]]$mses
    best_n_comps[[sim_group]][, n] <- results[[n]]$best_k
    train_conv_tracker[[sim_group]][[n]][["om_trate"]] <- results[[n]]$set_conv_train_track$om_trate
    train_conv_tracker[[sim_group]][[n]][["om_slog"]] <- results[[n]]$set_conv_train_track$om_slog
    train_conv_tracker[[sim_group]][[n]][["lcs"]] <- results[[n]]$set_conv_train_track$lcs
  }
}

method_names <- c("CFDA", "Windows", "OM-Trate (Hard)","OM-SLOG (Hard)", 
  "LCS (Hard)", "OM-Trate (Soft)","OM-SLOG (Soft)", "LCS (Soft)")

# Have to limit the data such that we make sure there is no issue with the 
# summarizing 
# should also report somewhere how many of each of the soft clustering are actually calculated in the average 

# 2. Transform the list and map the names
plot_data <- map_df(names(mse.sims), function(name) {
  matrix_data <- mse.sims[[name]]
  
  # Calculate row means - only works with n > 1
  row_means <- sqrt(rowMeans(matrix_data))
  
  # Create the data frame
  tibble(
    Method = factor(method_names, levels = method_names), # Keeps them in your specific order
    Mean_MSE = row_means,
    Simulation_Type = name
  )
})

# 3. Create the plot
plot_name <- paste0(sim_type, "plot.pdf")
  
pdf(plot_name,width=9,height=7)
ggplot(plot_data, aes(x = Method, y = Mean_MSE, color = Simulation_Type, group = Simulation_Type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) + # Makes text a bit more readable
  labs(
    title = "Mean RMSE per Method - Hard Difficulty",
    x = "Method",
    y = "Average RMSE (of Best Number of Components)",
    color = "Simulation Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Tilts labels if names are long
dev.off()

# Save data from full run 
saveRDS(mse.sims, paste0(sim_type,"_full_mses.rds"))
saveRDS(best_n_comps, paste0(sim_type,"_best_n_comps.rds"))
saveRDS(train_conv_tracker, paste0(sim_type, "_train_conv_tracker.rds"))
saveRDS(plot_data, paste0(sim_type, "_RMSE_performance.rds"))

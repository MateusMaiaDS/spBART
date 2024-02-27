# All BART-lite interaction
all_spbart_lite_interaction <- function(cv_element,
                                        nIknots_,
                                        node_min_size_,
                                        ntree_,
                                        seed_,
                                        j,
                                        alpha_,
                                        dif_order_,
                                        y_scale_,
                                        n_mcmc_,
                                        n_burn_,
                                        pen_basis_){

  # Doing a warming for the case whichI don't have
  # if(ntree_<50){
  #   stop("Use the all_bart() function instead.")
  # }

  if(length(cv_element)==7){
    kfold_bool <- TRUE
  } else {
    kfold_bool <- FALSE
  }
  # To replicate the results
  set.seed(seed_)

  if(kfold_bool){

    train <- cv_element$data_train
    test <- cv_element$data_test
    x_train <- cv_element$x_train
    x_test <- cv_element$x_test
    y_train <- cv_element$y_train
    y_test <- cv_element$y_test

  } else {
    train <- cv_element$train
    test <- cv_element$test

    # Getting the training elements
    x_train <- train[, colnames(train)!="y", drop = FALSE]
    x_test <- test[, colnames(train)!="y", drop = FALSE]
    y_train <- train %>% dplyr::pull("y")
    y_test <- test %>% dplyr::pull("y")
  }


  # Initialising df
  comparison_metrics <- data.frame(metric = NULL, value = NULL, model = NULL,fold = NULL)


  spBART_interaction <- rspBART(x_train = x_train,
                                x_test = x_test,y_train = y_train,
                                n_mcmc = n_mcmc_,node_min_size = node_min_size_,
                                alpha = alpha_,
                                n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
                                use_bs = FALSE,all_var = FALSE,
                                stump = FALSE,dif_order = dif_order_,scale_bool = y_scale_,
                                motrbart_bool = FALSE,
                                scale_init = FALSE,
                                interaction_term = TRUE,main_effects_pred = FALSE,
                                update_tau_beta = TRUE,
                                linero_sampler = FALSE,plot_preview = FALSE,
                                use_D_bool = FALSE,scale_basis_function = FALSE,
                                store_tree_fit = FALSE,varimportance_bool = TRUE,
                                robust_prior = FALSE,pen_basis = pen_basis_,eta = 1e-6,
                                center_basis = TRUE)


  n_burn_ <- n_burn_
  n_mcmc_ <- spBART_interaction$mcmc$n_mcmc

  # Calculating metrics for splinesBART
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = colMeans(spBART_interaction$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         train$y),
                                                            model = "psBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = colMeans(spBART_interaction$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         test$y),
                                                            model = "psBART",fold = j))

  # Calculating metrics for splinesBART
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(x = colMeans(spBART_interaction$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                        train$y),
                                                            model = "psBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(x = colMeans(spBART_interaction$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                        test$y),
                                                            model = "psBART",fold = j))


  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = colMeans(spBART_interaction$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         sds = rep(mean(spBART_interaction$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(train$y)))$CRPS,
                                                            model = "psBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = colMeans(spBART_interaction$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         sds = rep(mean(spBART_interaction$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(test$y)))$CRPS,
                                                            model = "psBART",fold = j))


  rm(spBART_interaction)



  return(comparison_metrics)

}



# Summarising all the metrics and results
wrapping_comparison <- function(result_){

  # Initialising df
  comparison_metrics <- data.frame(metric = NULL, value = NULL, model = NULL,fold = NULL)

  for(j in 1:length(result_)){


    n_burn_ <- result_[[j]]$spBART$mcmc$n_burn
    n_mcmc_ <- result_[[j]]$spBART$mcmc$n_mcmc

    if(!is.null(result_[[j]]$spBART_all)){
      # Calculating metrics for splinesBART
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(result_[[j]]$spBART_all$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             result_[[j]]$cv$train$y),
                                                                model = "spBART_all",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = colMeans(result_[[j]]$spBART_all$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             result_[[j]]$cv$test$y),
                                                                model = "spBART_all",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = result_[[j]]$cv$train$y ,
                                                                             means = colMeans(result_[[j]]$spBART_all$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             sds = rep(mean(result_[[j]]$spBART_all$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$train$y)))$CRPS,
                                                                model = "spBART_all",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = result_[[j]]$cv$test$y ,
                                                                             means = colMeans(result_[[j]]$spBART_all$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             sds = rep(mean(result_[[j]]$spBART_all$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$test$y)))$CRPS,
                                                                model = "spBART_all",fold = j))
    }

    # Calculating metrics for splinesBART
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(colMeans(result_[[j]]$spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           result_[[j]]$cv$train$y),
                                                              model = "spBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(colMeans(result_[[j]]$spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           result_[[j]]$cv$test$y),
                                                              model = "spBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = result_[[j]]$cv$train$y ,
                                                                           means = colMeans(result_[[j]]$spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           sds = rep(mean(result_[[j]]$spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$train$y)))$CRPS,
                                                              model = "spBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = result_[[j]]$cv$test$y ,
                                                                           means = colMeans(result_[[j]]$spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           sds = rep(mean(result_[[j]]$spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$test$y)))$CRPS,
                                                              model = "spBART",fold = j))

    # ============================
    # Calculating metrics for BART
    # ============================

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(result_[[j]]$bartmod$yhat.train.mean,
                                                                           result_[[j]]$cv$train$y),
                                                              model = "BART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(result_[[j]]$bartmod$yhat.test.mean,
                                                                           result_[[j]]$cv$test$y),
                                                              model = "BART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = result_[[j]]$cv$train$y ,
                                                                           means = result_[[j]]$bartmod$yhat.train.mean,
                                                                           sds = rep(mean(result_[[j]]$bartmod$sigma), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                              model = "BART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = result_[[j]]$cv$test$y ,
                                                                           means = result_[[j]]$bartmod$yhat.test.mean,
                                                                           sds = rep(mean(result_[[j]]$bartmod$sigma), length(result_[[j]]$cv$test$y) ))$CRPS,
                                                              model = "BART",fold = j))


    # ============================
    # Calculating metrics for softBART
    # ============================

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(result_[[j]]$softbartmod$y_hat_train_mean,
                                                                           result_[[j]]$cv$train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(result_[[j]]$softbartmod$y_hat_test_mean,
                                                                           result_[[j]]$cv$test$y),
                                                              model = "softBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = result_[[j]]$cv$train$y ,
                                                                           means = result_[[j]]$softbartmod$y_hat_train_mean,
                                                                           sds = rep(mean(result_[[j]]$softbartmod$sigma), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = result_[[j]]$cv$test$y ,
                                                                           means = result_[[j]]$softbart$y_hat_test_mean,
                                                                           sds = rep(mean(result_[[j]]$softbartmod$sigma), length(result_[[j]]$cv$test$y) ))$CRPS,
                                                              model = "softBART",fold = j))
    # ============================
    # Calculating metrics for MOTRBART
    # ============================


    if(result_[[j]]$spBART$prior$n_tree>1){
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                             result_[[j]]$cv$train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(colMeans(result_[[j]]$motrbart_pred),
                                                                             result_[[j]]$cv$test$y),
                                                                model = "motrBART",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = result_[[j]]$cv$train$y ,
                                                                             means = colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(result_[[j]]$motrbartmod$sigma2)), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = result_[[j]]$cv$test$y ,
                                                                             means = colMeans(result_[[j]]$motrbart_pred),
                                                                             sds = rep(mean(sqrt(result_[[j]]$motrbartmod$sigma2)), length(result_[[j]]$cv$test$y) ))$CRPS,
                                                                model = "motrBART",fold = j))
    } else {
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                             result_[[j]]$cv$train$y),
                                                                model = "motrBART",fold = j))


      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = result_[[j]]$cv$train$y ,
                                                                             means = colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(result_[[j]]$motrbartmod$sigma2)), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

    }

  }

  return(comparison_metrics)

}


# Getting a model to evaluate variable importance
var_importance_counter <- function(result_,rep_, only_sp_ = FALSE){

  if(!only_sp_){
    # Getting a counter for times that a variable is used in within a tree
    p_counter <- numeric(ncol(result_[[rep_]]$cv$train)-1)
    spBART <- result_[[rep_]]$spBART
  } else {
    spBART <- result_
    p_counter <- numeric(NCOL(spBART$data$x_train))
  }
  for(i in 501:spBART$mcmc$n_mcmc){

    for(t in 1:spBART$prior$n_tree){
      curr_tree <- spBART$mcmc$all_trees[[i]][[t]]
      terminals <- get_terminals(curr_tree)
      for(ell in 1:length(terminals)){
        p_counter[unique(curr_tree[[terminals[ell]]]$ancestors)] <- p_counter[unique(curr_tree[[terminals[ell]]]$ancestors)] + 1
      }
    }
  }

  return(round(p_counter/sum(p_counter),digits = 5))

}

# Getting a model to evaluate variable importance
tree_length_counter <- function(result_,rep_, only_sp = FALSE){

  # Getting a counter for times that a variable is used in within a tree
  if(!only_sp){
    spBART <- result_[[rep_]]$spBART
  } else {
    spBART <- result_
  }
  matrix_tree <- matrix(0,ncol = spBART$prior$n_tree, nrow = 2000)
  curr <- 0
  for(i in 501:spBART$mcmc$n_mcmc){
    curr <- curr + 1
    for(t in 1:spBART$prior$n_tree){
      curr_tree <- spBART$mcmc$all_trees[[i]][[t]]
      terminals <- get_terminals(curr_tree)
      matrix_tree[curr,t] <- length(terminals)
    }
  }

  return(matrix_tree)

}

# =======
# k_fold
# =======
kfold <- function(data_,
                  nfold_ = 10,
                  seed_){

  # Always set a seed
  if(is.null(seed_)){
    stop("Insert a valid seed.")
  }

  # Setting a seed
  set.seed(seed_)

  # Getting the cv
  cv_object_ <- vector("list",length = nfold_)

  # Getting the indentifier
  fold_index <- dismo::kfold(x = data_,k = nfold_)

  for(i in 1:nfold_){

    # Creating a list to store all information
    data_train <- data_[fold_index!=i,]
    data_test <- data_[fold_index==i,]
    x_train <- data_train[,which(colnames(data_train)!="y"), drop = FALSE]
    x_test <- data_test[,which(colnames(data_train)!="y"), drop = FALSE]
    y_train <- data_train[["y"]]
    y_test <- data_test[["y"]]

    # List witha all elements
    obj_ <- list(data_train = data_train,
                 data_test = data_test,
                 x_train = x_train,
                 x_test = x_test,
                 y_train = y_train,
                 y_test = y_test,
                 seed = seed_)

    cv_object_[[i]] <- obj_

  }

  # Return the list with all cv_objects
  return(cv_object_)

}

# spBART competitors
competitors_comparison_ <- function(cv_element,
                                    fold_,
                                    seed_,
                                    return_models = FALSE){

  # Setting a seed
  set.seed(seed_)

  # Initialising df
  comparison_metrics <- data.frame(metric = NULL,
                                   value = NULL,
                                   model = NULL,
                                   fold = NULL)


  if(length(cv_element)==7){
    kfold_bool <- TRUE
  } else {
    kfold_bool <- FALSE
  }
  # To replicate the results
  set.seed(seed_)

  if(kfold_bool){

    train <- cv_element$data_train
    test <- cv_element$data_test
    x_train <- cv_element$x_train
    x_test <- cv_element$x_test
    y_train <- cv_element$y_train
    y_test <- cv_element$y_test

  } else {
    train <- cv_element$train
    test <- cv_element$test

    # Getting the training elements
    x_train <- train[, colnames(train)!="y", drop = FALSE]
    x_test <- test[, colnames(train)!="y", drop = FALSE]
    y_train <- train %>% dplyr::pull("y")
    y_test <- test %>% dplyr::pull("y")
  }

  j <- fold_

  # Loading packages
  library(dbarts)
  library(SoftBart)
  library(earth)
  library(mgcv)
  library(MOTRbart)
  library(BASS)

  # Running all models
  bartmod <- dbarts::bart(x.train =  x_train,y.train = y_train,x.test = x_test)
  softbartmod <- SoftBart::softbart(X = x_train,Y = y_train,X_test = x_test)
  motrbartmod <- MOTRbart::motr_bart(x = x_train,y = y_train)
  motrbart_pred <- predict_motr_bart(motrbartmod,x_test,type = "mean")
  bassmod <- BASS::bass(xx = x_train,y = y_train)
  bass_pred <- predict(bassmod,x_test)

  marsmod <- earth(y~., data = train,degree = 2)
  marsmod_pred <- predict(marsmod,x_test)

  main_effects_formula <- paste0(paste0("s(",colnames(x_train),", bs = 'cr')"),collapse = "+")
  mgcv_formula <- formula(paste0("y~",paste0(c(main_effects_formula),collapse = "+")))

  gammod <- gam(mgcv_formula,data = train)
  gammod_pred <- predict(gammod,x_test)

  # ==================================
  # Calculating all metrics from now
  # ==================================

  # ============================
  # Calculating metrics for BART
  # ============================

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(bartmod$yhat.train.mean,
                                                                         train$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(bartmod$yhat.test.mean,
                                                                         test$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(bartmod$yhat.train.mean,
                                                                        train$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(bartmod$yhat.test.mean,
                                                                        test$y),
                                                            model = "BART",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = bartmod$yhat.train.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(train$y) ))$CRPS,
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = bartmod$yhat.test.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(test$y) ))$CRPS,
                                                            model = "BART",fold = j))


  # ============================
  # Calculating metrics for softBART
  # ============================

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(softbartmod$y_hat_train_mean,
                                                                         train$y),
                                                            model = "softBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(softbartmod$y_hat_test_mean,
                                                                         test$y),
                                                            model = "softBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(softbartmod$y_hat_train_mean,
                                                                        train$y),
                                                            model = "softBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(softbartmod$y_hat_test_mean,
                                                                        test$y),
                                                            model = "softBART",fold = j))


  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = softbartmod$y_hat_train_mean,
                                                                         sds = rep(mean(softbartmod$sigma), length(train$y) ))$CRPS,
                                                            model = "softBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = softbartmod$y_hat_test_mean,
                                                                         sds = rep(mean(softbartmod$sigma), length(test$y) ))$CRPS,
                                                            model = "softBART",fold = j))
  # ============================
  # Calculating metrics for MOTRBART
  # ============================
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(colMeans(motrbartmod$y_hat),
                                                                         train$y),
                                                            model = "motrBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(motrbart_pred,
                                                                         test$y),
                                                            model = "motrBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(colMeans(motrbartmod$y_hat),
                                                                        train$y),
                                                            model = "motrBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(motrbart_pred,
                                                                        test$y),
                                                            model = "motrBART",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = colMeans(motrbartmod$y_hat),
                                                                         sds = rep(mean(sqrt(motrbartmod$sigma2)), length(train$y) ))$CRPS,
                                                            model = "motrBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = motrbart_pred,
                                                                         sds = rep(mean(sqrt(motrbartmod$sigma2)), length(test$y) ))$CRPS,
                                                            model = "motrBART",fold = j))

  # ============================
  # Calculating metrics for MARS
  # ============================
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(c(marsmod$fitted.values),
                                                                         train$y),
                                                            model = "MARS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(c(marsmod_pred),
                                                                         test$y),
                                                            model = "MARS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(c(marsmod$fitted.values),
                                                                        train$y),
                                                            model = "MARS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(c(marsmod_pred),
                                                                        test$y),
                                                            model = "MARS",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = c(marsmod$fitted.values),
                                                                         sds = rep(sigma(marsmod), length(train$y) ))$CRPS,
                                                            model = "MARS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = c(marsmod_pred),
                                                                         sds = rep(sigma(marsmod), length(test$y) ))$CRPS,
                                                            model = "MARS",fold = j))

  # ============================
  # Calculating metrics for GAM
  # ============================
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(c(gammod$fitted.values),
                                                                         train$y),
                                                            model = "GAM",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(c(gammod_pred),
                                                                         test$y),
                                                            model = "GAM",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(c(gammod$fitted.values),
                                                                        train$y),
                                                            model = "GAM",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(c(gammod_pred),
                                                                        test$y),
                                                            model = "GAM",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = c(gammod$fitted.values),
                                                                         sds = rep(sigma(marsmod), length(train$y) ))$CRPS,
                                                            model = "GAM",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = c(gammod_pred),
                                                                         sds = rep(sigma(marsmod), length(test$y) ))$CRPS,
                                                            model = "GAM",fold = j))


  # =======================================
  #     Doing a comparison with BASS
  # =======================================

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(c(bassmod$yhat.mean),
                                                                         train$y),
                                                            model = "BASS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(colMeans(bass_pred),
                                                                         test$y),
                                                            model = "BASS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(c(bassmod$yhat.mean),
                                                                        train$y),
                                                            model = "BASS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(colMeans(bass_pred),
                                                                        test$y),
                                                            model = "BASS",fold = j))


  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = c(bassmod$yhat.mean),
                                                                         sds = rep(mean(bassmod$s2), length(train$y) ))$CRPS,
                                                            model = "BASS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = colMeans(bass_pred),
                                                                         sds = rep(mean(bassmod$s2), length(test$y) ))$CRPS,
                                                            model = "BASS",fold = j))

  # Returning a list with all elements
  if(return_models){
    return(list(comparison_metrics = comparison_metrics,
                bartmod = bartmod,
                softbartmod = softbartmod,
                motrbartmod = motrbartmod,
                motrbart_pred = motrbart_pred,
                marsmod = marsmod,
                marsmod_pred = marsmod_pred,
                gammod = gammod,
                gammod_pred = gammod_pred,
                bassmod = bassmod,
                bass_pred = bass_pred))
  } else {
    return(list(comparison_metrics = comparison_metrics))
  }

}

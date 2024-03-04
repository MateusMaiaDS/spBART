# Cleaning the working directory space
rm(list=ls())
# Plotting a raw boxplot
library(ggplot2)
library(tidyverse)
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
library(doParallel)
source("R/sim_functions_27.R")
source("R/main_function_27.R")
source("R/cv_functions.R")


# Getting a single plot
# single_comparison_df <- wrapping_comparison(result_ = result)

# For the lite object

ntree_ <- 3
type_ <- "airquality"
# type_ <- "friedman_break"
# type_ <- "friedman"
# type_ <- "non.and.smooth.main.formula"

nIknots_ <- 20
node_min_size_ <- 25
n_mcmc_ <- 10000
n_burn_ <- 5000
dif_order_ <- 2

# Retrieving rpsBART results
# Getting the for multiple knots
rep_vec <- 1:10
seed_ <- 43
all_boxplot_df <- data.frame(metric = NULL,
                             value = NULL,
                             model = NULL,
                             fold = NULL,
                             n = NULL)
for(jj in rep_vec){
  rps_bart_result <- readRDS(paste0("~/spline_bart_lab/preliminar_results/rspBART27/airquality/v31_rep_",
                                    jj,"_intercept_psBART_seed_",seed_,"_nIknots_",nIknots_,
                                    "_ntree_",ntree_,"_alpha_0.5_dif_2_nmin_",node_min_size_,"_nmcmc_",n_mcmc_,"_nburn_",n_burn_,".Rds"))
  all_boxplot_df <- rbind(all_boxplot_df,rps_bart_result)

}



# Retrieving the results from its cmpetitors
competitors_result <- readRDS(paste0("~/spline_bart_lab/preliminar_results/rspBART27/",type_,"/competitors_seed_",seed_,".Rds"))
competitors_df <- do.call(rbind,competitors_result) %>% do.call(rbind,.)

all_boxplot_df <- rbind(all_boxplot_df,competitors_df)

all_boxplot_df <- all_boxplot_df %>%
  dplyr::mutate(model = factor(model, c("psBART","softBART","MARS","BASS","GAM","BART","motrBART"))) %>%
  dplyr::mutate(metric = factor(metric, c("rmse_train","mae_train","crps_train","rmse_test","mae_test","crps_test")))
fff <- 4
par(mfrow = c(1,1))

# ntree_ <- result[[fff]]$
ranking <- all_boxplot_df %>% filter(metric=="rmse_test") %>% group_by(fold) %>% mutate(rank = rank(value)) %>%
  ungroup(fold) %>% group_by(model) %>% summarise(mean_rank = mean(rank))

ranking %>% arrange(mean_rank)


# single_comparison_df <- do.call(rbind,result)
ggplot(all_boxplot_df)+
  facet_wrap(~metric, scales = "free_y")+
  geom_boxplot(mapping = aes(x = model, y = value))+
  ggtitle(paste0("Simulation Scenario: ",type_))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

all_boxplot_df %>% filter(metric=="rmse_test") %>% group_by(model) %>%
  summarise(mean_rmse = mean(value)) %>% arrange(mean_rmse)

# result %>% do.call(rbind,.) %>% filter(model=="spBART_inter") %>% group_by(metric) %>% summarise(mean_val = mean(value))
# result2 %>% do.call(rbind,.) %>% filter(model=="spBART_all") %>% group_by(metric) %>% summarise(mean_val = mean(value))

# result[[fff]]$cv$train %>% plot
# points(result[[fff]]$spBART$data$x_train$x,result[[fff]]$spBART$y_train_hat %>% colMeans(),pch=20)
# result[[fff]]$spBART$all_tau %>% plot(type = "l")
# points(result[[1]]$spBART$data$x_train$x,result[[1]]$softbartmod$y_hat_train_mean ,pch=20, col = "blue")
# points(result[[fff]]$spBART$data$x_train$x,result[[fff]]$motrbartmod$y_hat %>% colMeans() ,pch=20, col = "orange")

par(mfrow= c(1,3))
plot(result[[fff]]$cv$train$x,result[[fff]]$cv$train$y, main = paste0("spBART - Ntree = ",ntree_ ), xlab = "x", ylab = "y")
points(result[[fff]]$spBART$data$x_train$x,result[[fff]]$spBART$y_train_hat %>% colMeans(),pch=20)

# If we are interested into the quantiles
# quantiles_ <- result[[fff]]$spBART$y_train_hat %>% apply(2,function(x){quantile(x,probs = c(0.025,0.975))})
# points(result[[fff]]$spBART$data$x_train$x,quantiles_[1,] ,pch=20,col = ggplot2::alpha("red",0.1))
# points(result[[fff]]$spBART$data$x_train$x,quantiles_[2,] ,pch=20,col = ggplot2::alpha("red",0.1))

plot(result[[fff]]$cv$train$x,result[[fff]]$cv$train$y, main = paste0("softBART - Ntree = ",20), xlab = "x", ylab = "y")
points(result[[fff]]$spBART$data$x_train$x,result[[fff]]$softbartmod$y_hat_train_mean ,pch=20, col = "blue")

# Case for all variables comparison
# plot(result[[fff]]$cv$train$x,result[[fff]]$cv$train$y, main = paste0("spBART_all - Ntree = ",20), xlab = "x", ylab = "y")
# points(result[[1]]$spBART_all$data$x_train$x,result[[1]]$spBART_all$y_train_hat %>% colMeans() ,pch=20, col = "darkgreen")

plot(result[[fff]]$cv$train$x,result[[fff]]$cv$train$y, main = paste0("MOTRBART - Ntree = ",ntree_ ), xlab = "x", ylab = "y")
points(result[[fff]]$spBART$data$x_train$x,result[[fff]]$motrbartmod$y_hat %>% colMeans() ,pch=20, col = "orange")


# ==== PLOTTING FOR THE TEST SET ==================

par(mfrow= c(1,3))
plot(result[[fff]]$cv$test$x,result[[fff]]$cv$test$y, main = paste0("spBART - Ntree = ",ntree_ ), xlab = "x", ylab = "y")
points(result[[fff]]$spBART$data$x_test$x,result[[fff]]$spBART$y_test_hat %>% colMeans(),pch=20)

plot(result[[fff]]$cv$test$x,result[[fff]]$cv$test$y, main = paste0("softBART - Ntree = ",20), xlab = "x", ylab = "y")
points(result[[fff]]$spBART$data$x_test$x,result[[fff]]$softbartmod$y_hat_test_mean ,pch=20, col = "blue")

# Case for all variables comparison
# plot(result[[fff]]$cv$test$x,result[[fff]]$cv$test$y, main = paste0("spBART_all - Ntree = ",20), xlab = "x", ylab = "y")
# points(result[[1]]$spBART_all$data$x_test$x,result[[1]]$spBART_all$y_test_hat %>% colMeans() ,pch=20, col = "darkgreen")

plot(result[[fff]]$cv$test$x,result[[fff]]$cv$test$y, main = paste0("MOTRBART - Ntree = ",ntree_ ), xlab = "x", ylab = "y")
points(result[[fff]]$spBART$data$x_test$x,result[[fff]]$motrbartmod$y_hat %>% colMeans() ,pch=20, col = "orange")



# Plotting all results summary
vec_nIknots <- c(2,5)
vec_ntree <- c(1,10,20,50)
vec_alpha <- c(0.95)
allvar_ <- TRUE
update_beta_ <- TRUE
dif_order_ <- 1
inter_ <- TRUE
# Selecting a simulated scenarion
# (1): "oned_break" one dimensionnal sin(2*x) with a break
# (2): "friedman_nointer_nonoise": four-dimensional friedmna setting with no interaction terms and no extra X noise variables
type_ <- c("friedman_inter_noise")

all_comparison_df <- data.frame(metric = NULL,
                                value = NULL,
                                model = NULL,
                                fold = NULL,
                                nIknots = NULL,
                                ntree = NULL,
                                alpha = NULL)

for(m in 1:length(vec_alpha)){
  for(k in 1:length(vec_ntree)){
    for(j in 1:length(vec_nIknots)){

      if(type_ == "friedman_nointer_nonoise"){
        result <- readRDS(paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART8/friedman_noint_nonoise/oned_n_250_sd_1_nIknots_",vec_nIknots[j],"_ntree_",vec_ntree[k],"_bs_FALSE_motr_bart_FALSE_allvar_",allvar_,"_stump_TRUE_sinit_FALSE_alpha_",vec_alpha[m],"_uptaubeta_TRUE.Rds"))
      } else if(type_ == "oned_break"){
        result <- readRDS(paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART8/oned_n_250_sd_1_nIknots_",vec_nIknots[j],"_ntree_",vec_ntree[k],"_bs_FALSE_motr_bart_FALSE_alpha_",vec_alpha[m],".Rds"))
      } else if(type_ == "friedman_nointer_noise"){
        result <- readRDS(paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART8/friedman_noint_noise/oned_n_250_sd_1_nIknots_",vec_nIknots[j],"_ntree_",vec_ntree[k],"_bs_FALSE_motr_bart_FALSE_allvar_",allvar_,"_stump_TRUE_sinit_FALSE_alpha_0.95_uptaubeta_TRUE_dif_1.Rds"))
      } else if(type_ == "friedman_inter_noise"){
        if(vec_ntree[k]==50 & vec_nIknots[j]==10){
          next
        }
        result <- readRDS(paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART8/friedman/new_interaction_",inter_,"_oned_n_250_sd_1_nIknots_",vec_nIknots[j],"_ntree_",vec_ntree[k],"_bs_FALSE_motr_bart_FALSE_allvar_",allvar_,"_stump_TRUE_sinit_TRUE_alpha_0.95_uptaubeta_",update_beta_,"_dif_",dif_order_,".Rds"))
      }
      all_comparison_df <- result %>% do.call(rbind,.) %>% mutate(nIknots = vec_nIknots[j], ntree = vec_ntree[k],alpha = vec_alpha[m]) %>% rbind(all_comparison_df)
      print(paste0("Running iteration .... : ", j))
      rm(result)

    }
    print(paste0("Running iteration tree .... : ", k))

  }
  print(paste0("Running alpha .... : ", m))
}

# PLotting results

n_rep_ <- 10
alpha_ <- 0.95

all_comparison_df  %>% filter( alpha == alpha_ ) %>%
  dplyr::mutate(nIknots = as.factor(nIknots), ntree = as.factor(ntree),alpha = as.factor(alpha)) %>%
  filter((model == "motrBART" & nIknots == "2" ) |
           (model == "BART" & nIknots == "2" & ntree == "10" ) |
           (model == "softBART" & nIknots == "2" & ntree == "10") |
           (model == "spBART_inter") |
           (model == "spBART_all"))   %>%
  filter(metric == "rmse_train") %>%
  ggplot()+
  ggtitle(paste0(n_rep_,"-fold CV // Alpha (spBART) = ",alpha_," // Facet:(n_tree)"),subtitle = type_)+
  ylab("RMSE train set") +
  facet_wrap(~ntree)+
  geom_boxplot(mapping = aes(x = model, y = value, col = nIknots))+
  theme_bw()

all_comparison_df %>%  ggplot()+
  facet_wrap(~metric)+
  geom_boxplot(mapping = aes(x = model, y = value))

# Keefe experiment for different SD's
sd_vec <- c(1,5,10,25)
n_tree_vec <- c(10,50)

keefe_comparison <- data.frame(metric = NULL,
                               value = NULL,
                               model = NULL,
                               fold = NULL,
                               ntree = NULL,
                               sd = NULL)

for(i in 1:length(sd_vec)){
  for(j in 1:length(n_tree_vec)){
    result <- readRDS(paste0("~/spline_bart_lab/preliminar_results/rspBART12/friedman/v3_new_interaction_TRUE_oned_n_250_sd_",sd_vec[i],"_nIknots_2_ntree_",n_tree_vec[j],"_bs_FALSE_motr_bart_FALSE_allvar_TRUE_stump_FALSE_sinit_FALSE_alpha_0.5_uptaubeta_TRUE_dif_1.Rds"))
    keefe_comparison<- result %>% do.call(rbind,.) %>% mutate(ntree = n_tree_vec[j],sd = sd_vec[i]) %>% rbind(keefe_comparison)
    print(paste0("Running iteration .... : ", j))
  }

}


keefe_comparison %>%
  dplyr::mutate(ntree = as.factor(ntree),sd = as.factor(sd)) %>%
  filter((model == "motrBART"  ) |
           (model == "BART"  & ntree == "10" ) |
           (model == "softBART" & ntree == "10") |
           (model == "spBART_inter")) %>%  # |
  # (model == "spBART_all"))   %>%
  filter(metric == "crps_test") %>%
  ggplot()+
  ggtitle(paste0(n_rep_,"-fold CV // Alpha (spBART) = ",alpha_," // Facet:(sd)"),subtitle = type_)+
  ylab("CRPS test set") +
  facet_wrap(~sd, scales = "free_y")+
  geom_boxplot(mapping = aes(x = model, y = value, col = ntree))+
  theme_bw()


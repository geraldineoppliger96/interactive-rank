sudo_wilcoxon = function(dat, n_permute, alg_type, rank_type, max_nodes = NULL, paired = FALSE){
  if (alg_type == "sudo") {
    stat_func = function(x) {sudo_RF(x, rank_type = rank_type, max_nodes = max_nodes)}
  } else if (alg_type == "mask") {
    stat_func = function(x) {mask_RF(x, rank_type = rank_type, max_nodes = max_nodes)}
  } else if (alg_type == "OB") {
    stat_func = function(x) {OB_RF(x, rank_type = rank_type, abs_type = FALSE)}
  } else if (alg_type == "semi") {
    stat_func = function(x) {semi_RF(x, rank_type = rank_type, abs_type = FALSE)}
  } else if (alg_type == "OB-abs") {
    stat_func = function(x) {OB_RF(x, rank_type = rank_type, abs_type = TRUE)}
  } else if (alg_type == "semi-abs") {
    stat_func = function(x) {semi_RF(x, rank_type = rank_type, abs_type = TRUE)}
  }
  
  W_obs = stat_func(dat); print(W_obs)
  n = length(dat$Z)
  W_permute = foreach(i = 1:n_permute, .combine = cbind,
                      .export = c("sudo_RF", "mask_RF", "OB_RF", "semi_RF"),
                      .packages = c("foreach", "randomForest")) %dopar% {
    dat_permute = dat; 
    if (paired) {
      temp = base::sample(dat$Z[1:(n/2)])
      dat_permute$Z = c(temp, 1 - temp)
    } else {
      dat_permute$Z = base::sample(dat_permute$Z)
    }
    stat_func(dat_permute)
                      }
  # df = data.frame(x = as.vector(W_permute))
  # p = ggplot(df, aes(x=x)) +
  #   geom_histogram(aes(y=..density..), colour="black", fill="white")+
  #   geom_density(alpha=.2, fill="#FF6666") +
  #   xlab("permuted values: pseudo-Wilcoxon-value") +
  #   geom_vline(aes(xintercept=W_obs),
  #              color="blue", linetype="dashed", size=1)
  # plot(p)
  # ggsave(filename = paste("figure/permute_oob_alter.png", sep = ""), plot = p,
  #        width = 4, height = 3.6)
  pval = mean(W_permute > W_obs)
  return(pval)
}

sudo_RF = function(dat, rank_type, max_nodes) {
  n = nrow(dat); ncov = ncol(dat) - 1
  model = randomForest(Y ~ ., data = dat, maxnodes = max_nodes)
  yhat_true = predict(model)
  dat_false = dat; dat_false$Z = 1 - dat_false$Z
  yhat_false = predict(model, newdata = dat_false)
  
  E = abs(yhat_false - dat$Y) - abs(yhat_true - dat$Y)
  if (rank_type) {
    W = sum((2*(E > 0) - 1)*rank(abs(E)))
  } else {
    W = sum(E)
  }
  return(W)
}

mask_RF = function(dat, rank_type, max_nodes) {
  n = nrow(dat); ncov = ncol(dat) - 1
  model_full = randomForest(Y ~ ., data = dat, maxnodes = max_nodes)
  yhat_full = predict(model_full)
  model_mask = randomForest(Y ~ .-Z, data = dat, maxnodes = max_nodes)
  yhat_mask = predict(model_mask)
  
  E = abs(yhat_mask - dat$Y) - abs(yhat_full - dat$Y) 
  if (rank_type) {
    W = sum((2*(E > 0) - 1)*rank(abs(E)))
  } else {
    W = sum(E)
  }
  return(W)
}

OB_RF = function(dat, rank_type, abs_type) {
  model = randomForest(Y ~ ., data = dat)
  dat_false = dat; dat_false$Z = 1 - dat_false$Z
  yhat_false = predict(model, newdata = dat_false)
  if (rank_type) {
    if (abs_type) {
      E = (dat$Y - yhat_false)*(2*dat$Z - 1)
      W = sum((2*(E > 0) - 1)*rank(abs(E)))
    } else {
      W = sum(rank(dat$Y - yhat_false)*(2*dat$Z - 1))
    }
  } else {
    if (abs_type) {
      W = sum(abs(dat$Y - yhat_false))
    } else {
      W = sum((dat$Y - yhat_false)*(2*dat$Z - 1))
    }
  }
  return(W)
}

semi_RF = function(dat, rank_type, abs_type) {
  model = randomForest(Y ~ ., data = dat)
  yhat_true = predict(model)
  dat_false = dat; dat_false$Z = 1 - dat_false$Z
  yhat_false = predict(model, newdata = dat_false)
  if (rank_type) {
    if (abs_type) {
      E = (2*dat$Y - yhat_false - yhat_true)*(2*dat$Z - 1)
      W = sum((2*(E > 0) - 1)*rank(abs(E)))
    } else {
      W = sum(rank(2*dat$Y - yhat_false - yhat_true)*(2*dat$Z - 1))
    }
  } else {
    if (abs_type) {
      W = sum(abs(2*dat$Y - yhat_false - yhat_true))
    } else {
      W = sum((2*dat$Y - yhat_false - yhat_true)*(2*dat$Z - 1))
    }
  }
  return(W)
}



permute_wilcoxon = function(Y_diff, Z_diff, n_permute){
  stat_func = function(Y_diff, Z_diff){
    sum(Y_diff*Z_diff)
  }
  W_obs = stat_func(Y_diff, Z_diff); print(W_obs)
  
  W_permute = foreach(i = 1:n_permute, .combine = cbind) %dopar% {
                        Z_diff_permute = base::sample(Z_diff)
                        stat_func(Y_diff, Z_diff_permute)
  }
  pval = mean(W_permute > W_obs)
  return(pval)
}


sudo_wilcoxon_mixed = function(dat, n_pair, n_permute, alg_type, rank_type, max_nodes = NULL){
  if (alg_type == "sudo") {
    stat_func = function(x) {sudo_RF(x, rank_type = rank_type, max_nodes = max_nodes)}
  } else if (alg_type == "mask") {
    stat_func = function(x) {mask_RF(x, rank_type = rank_type, max_nodes = max_nodes)}
  }
  
  W_obs = stat_func(dat); print(W_obs); n = length(dat$Z)
  W_permute = foreach(i = 1:n_permute, .combine = cbind,
                      .export = c("sudo_RF", "mask_RF"),
                      .packages = c("foreach", "randomForest")) %dopar% {
      dat_permute = dat;
      temp = base::sample(dat$Z[1:n_pair])
      dat_permute$Z = c(temp, 1 - temp, base::sample(dat$Z[(2*n_pair + 1):n]))
      stat_func(dat_permute)
    }
  pval = mean(W_permute > W_obs)
  return(pval)
}

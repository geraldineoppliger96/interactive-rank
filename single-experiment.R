experiment_pair_noside = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  rejection = matrix(ncol = 2, nrow = R); colnames(rejection) = c("wilcoxon", "martingale")
  for (i in 1:R) {
    print(i)
    Y = sample_pair_noside(n = n, p = p, delta = delta, D = D)
    if (all(Y$Y1 == Y$Y2)) {
      rejection[i, "wilcoxon"] = FALSE 
    } else {
      rejection[i, "wilcoxon"] = 
        wilcox.test(Y$Y1, Y$Y2, paired = TRUE)$p.value < alpha
    }
    rejection[i, "martingale"] = i_sign_rank(Y1 = Y$Y1, Y2 = Y$Y2, ub = ub, alpha = alpha)
  }
  return(rejection)
}

experiment_pair_side = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i) {
    print(i)
    rejection = vector(length = length(methods_pair)); names(rejection) = methods_pair
    
    dat = sample_pair_side(t = t, n = n, m = m, delta = delta, prob = prob_assign,
                           C_delta = C_delta, e_D = e_D, dat_type = dat_type)
    dat_merge = data.frame(Z = c(dat$Z1, dat$Z2), Y = c(dat$Y1, dat$Y2), X = rbind(dat$X1, dat$X2))
    dat_diff = data.frame(Z = (dat$Z1 - dat$Z2 + 1)/2, Y = dat$Y1 - dat$Y2, X = cbind(dat$X1, dat$X2))
    if (alg_type == "linear") {
      df = data.frame(Y = dat$Y1 - dat$Y2, X = dat$X1 - dat$X2)
      e = lm(Y ~ .^2 -1, data = df)$residuals
    } else if (alg_type == "nonlinear") {
      df = data.frame(Y = dat$Y1, X = dat$X1)
      e1 = df$Y - randomForest(Y ~., data = df)$predicted
      df = data.frame(Y = dat$Y2, X = dat$X2)
      e2 = df$Y - randomForest(Y ~., data = df)$predicted
      e = e1 - e2
    }
    
    if ("Wilcoxon" %in% methods_pair) {
      rejection["Wilcoxon"] = wilcox.test(dat$Y1, dat$Y2, paired = TRUE)$p.value
    }
    if ("Wilcoxon-value" %in% methods_pair) {
      rejection["Wil-value"] = 
        permute_wilcoxon(Y_diff = dat$Y1-dat$Y2, Z_diff = dat$Z1 - dat$Z2, n_permute = n_permute) 
    }
    if ("MannW" %in% methods_pair) {
      rejection["MannW"] = 
        wilcox.test(dat_diff$Y ~ dat_diff$Z, paired = FALSE)$p.value 
    }
    if ("covadj-Wilcoxon" %in% methods_pair) {
      rejection["covadj-Wilcoxon"] = 
        wilcox.test((dat$Z1 - dat$Z2)*e, rep(0, n), paired = TRUE)$p.value 
    }
    if ("covadj-Wilcoxon-value" %in% methods_pair) {
      rejection["covadj-Wilcoxon-value"] = 
        permute_wilcoxon(Y_diff = e, Z_diff = dat$Z1 - dat$Z2, n_permute = n_permute)
    }
    if ("covadj-MannW" %in% methods_pair) {
      rejection["covadj-MannW"] = 
        wilcox.test(e ~ dat_diff$Z, paired = FALSE)$p.value 
    }
    if ("indep-test" %in% methods_pair) {
      rejection["indep-test"] = dcov.test(e, dat$Z1 - dat$Z2, R = 100)$p.value
    }
    if ("cor-test" %in% methods_pair) {
      rejection["cor-test"] = cor.test(e, dat$Z1 - dat$Z2)$p.value
    }
    if ("martingale-pair-combine" %in% methods_pair) {
      rejection["martingale-pair-combine"] =
        split_wilcoxon(dat = dat, ub = ub, pred_func = combine_predict,
                       prob_assign = prob_assign)
    }
    if ("martingale-pair-sep" %in% methods_pair) {
      rejection["martingale-pair-sep"] =
        split_wilcoxon(dat = dat, ub = ub, pred_func = sep_predict,
                       prob_assign = prob_assign)
    }
    if ("pseudo-Wilcoxon-sep-rank" %in% methods_pair) {
      rejection["pseudo-Wilcoxon-sep-rank"] = 
        sudo_wilcoxon(dat_merge, n_permute = n_permute, paired = TRUE,
                      alg_type = "sudo", rank_type = TRUE, max_nodes = NULL) 
    }
    if ("pseudo-Wilcoxon-sep-value" %in% methods_pair) {
      rejection["pseudo-Wilcoxon-sep-value"] = 
        sudo_wilcoxon(dat_merge, n_permute = n_permute, paired = TRUE,
                      alg_type = "sudo", rank_type = FALSE, max_nodes = NULL) 
    }
    if ("pseudo-Wilcoxon-combine-rank" %in% methods_pair) {
      rejection["pseudo-Wilcoxon-combine-rank"] =
        sudo_wilcoxon(dat_diff, n_permute = n_permute,
                      alg_type = "sudo", rank_type = TRUE, max_nodes = NULL) 
    }
    if ("pseudo-Wilcoxon-combine-value" %in% methods_pair) {
      rejection["pseudo-Wilcoxon-combine-value"] =
        sudo_wilcoxon(dat_diff, n_permute = n_permute,
                      alg_type = "sudo", rank_type = FALSE, max_nodes = NULL)
    }
    if ("mask-Wilcoxon-sep" %in% methods_pair) {
      rejection["mask-Wilcoxon-sep"] = sudo_wilcoxon(dat_merge, n_permute = n_permute,
                                                    alg_type = "mask", max_nodes = NULL) 
    }
    if ("mask-Wilcoxon-combine" %in% methods_pair) {
      rejection["mask-Wilcoxon-combine"] = sudo_wilcoxon(dat_diff, n_permute = n_permute,
                                                    alg_type = "mask", max_nodes = NULL) 
    }
    if ("interactive" %in% methods_pair) {
      if(prob_assign == 1/2) {
        rejection["interactive"] = i_wilcoxon(X1 = dat$X1, X2 = dat$X2, Y1 = dat$Y1, Y2 = dat$Y2,
                                              Z1 = dat$Z1, Z2 = dat$Z2, ub = ub, alpha = alpha,
                                              dat_type = "linear", alg_type = "em-lm")
      } else {
        rejection["interactive"] =
          i_wilcoxon_two(X1 = dat$X1, X2 = dat$X2, Y1 = dat$Y1, Y2 = dat$Y2,
                         Z1 = dat$Z1, Z2 = dat$Z2, ub = ub,
                         alpha = alpha, prob_assign = prob_assign)
      }
      
    }
    return(rejection)
  }
  #rejections = mclapply(1:R, wrapper_func, mc.cores = 4)
  rejections = lapply(1:R, wrapper_func)
  return(rejections)
}



experiment_unpair = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    rejection = vector(length = length(methods_unpair)); names(rejection) = methods_unpair
    
    dat = sample_unpair(t = t, n = n, m = m, delta = delta, prob = prob_assign,
                        C_delta = C_delta, Cf = Cf, e_D = e_D, dat_type = dat_type,
                        redundant = redundant, blind = blind)
    df = data.frame(Y = dat$Y, dat$X);
    if (alg_type %in% c("linear")) {
      e = rlm(Y ~ .^2 + I(X3^2) - 1, data = df)$residuals
    } else if(alg_type == "RF"){
      e = dat$Y - randomForest(Y ~ ., data = df)$predicted
    } 
    
    if ("MannW" %in% methods_unpair) {
      rejection["MannW"] = wilcox.test(dat$Y ~ dat$Z, paired = FALSE)$p.value 
    }
    if ("Wil" %in% methods_unpair) {
      rejection["Wil"] =
        wilcox.test(dat$Y*(2*dat$Z - 1), rep(0, n), paired = TRUE)$p.value 
    }
    if ("Wil-value" %in% methods_unpair) {
      rejection["Wil-value"] = 
        permute_wilcoxon(Y_diff = dat$Y, Z_diff = (2*dat$Z - 1), n_permute = n_permute) 
    }
    if ("covadj-MannW" %in% methods_unpair) {
      rejection["covadj-MannW"] = wilcox.test(e ~ dat$Z, paired = FALSE)$p.value 
    }
    if ("covadj-Wil" %in% methods_unpair) {
      rejection["covadj-Wil"] = 
        wilcox.test((2*dat$Z - 1)*e, rep(0, n), paired = TRUE)$p.value 
    }
    if ("covadj-Wil-value" %in% methods_unpair) {
      rejection["covadj-Wil-value"] = 
        permute_wilcoxon(Y_diff = e, Z_diff = (2*dat$Z - 1), n_permute = n_permute) 
    }
    if ("CATE-test" %in% methods_unpair) {
      rejection["CATE-test"] = CATE_test_pval(X = dat$X, Z = dat$Z, Y = dat$Y)
    }
    if ("indep-test" %in% methods_unpair) {
      rejection["indep-test"] = dcov.test(e, dat$Z, R = 100)$p.value 
    }
    if ("cor-test" %in% methods_unpair) {
      rejection["cor-test"] = cor.test(e, dat$Z)$p.value
    }
    if ("split-MannW" %in% methods_unpair) {
      rejection["split-MannW"] =
        split_mann(X = dat$X, Z = dat$Z, Y = dat$Y, ub = ub, prob_assign = prob_assign) 
    }
    if ("interactive" %in% methods_unpair) {
      rejection["interactive"] = 
        i_mann_pval(X = dat$X, Z = dat$Z, Y = dat$Y, prob_assign = prob_assign)
    }
    if ("interactive-default" %in% methods_unpair) {
      rejection["interactive-default"] = 
        i_mann_pval(X = dat$X, Z = dat$Z, Y = dat$Y, prob_assign = prob_assign, type = "default")
    }
    if ("interactive-naive" %in% methods_unpair) {
      rejection["interactive-naive"] = 
        i_mann_pval(X = dat$X, Z = dat$Z, Y = dat$Y, prob_assign = prob_assign, type = "naive")
    }
    if ("order-martingale" %in% methods_unpair) {
      rejection["order-martingale"] = order_mann_pval(Z = dat$Z, Y = e)
    }
    if(0){
    if ("pseudo-Wil-rank" %in% methods_unpair) {
      rejection["pseudo-Wil-rank"] =
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute, 
                      rank_type = TRUE, alg_type = "sudo", max_nodes = NULL) 
    }
    if ("pseudo-Wil-value" %in% methods_unpair) {
      rejection["pseudo-Wil-value"] =
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute, 
                      rank_type = FALSE, alg_type = "sudo", max_nodes = NULL) 
    }
    if ("mask-Wil-rank" %in% methods_unpair) {
      rejection["mask-Wil-rank"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = TRUE, alg_type = "mask", max_nodes = NULL) 
    }
    if ("mask-Wil-value" %in% methods_unpair) {
      rejection["mask-Wil-value"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = FALSE, alg_type = "mask", max_nodes = NULL) 
    }
    if ("false-Wil-rank" %in% methods_unpair) {
      rejection["false-Wil-rank"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = TRUE, alg_type = "OB", max_nodes = NULL)
    }
    if ("false-Wil-value" %in% methods_unpair) {
      rejection["false-Wil-value"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = FALSE, alg_type = "OB", max_nodes = NULL)
    }
    if ("aver-Wil-rank" %in% methods_unpair) {
      rejection["aver-Wil-rank"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = TRUE, alg_type = "semi", max_nodes = NULL)
    }
    if ("aver-Wil-value" %in% methods_unpair) {
      rejection["aver-Wil-value"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = FALSE, alg_type = "semi", max_nodes = NULL)
    }
    if ("false-abs-Wil-rank" %in% methods_unpair) {
      rejection["false-abs-Wil-rank"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute, 
                      rank_type = TRUE, alg_type = "OB-abs", max_nodes = NULL)
    }
    if ("false-abs-Wil-value" %in% methods_unpair) {
      rejection["false-abs-Wil-value"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = FALSE, alg_type = "OB-abs", max_nodes = NULL)
    }
    if ("aver-abs-Wil-rank" %in% methods_unpair) {
      rejection["aver-abs-Wil-rank"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = TRUE, alg_type = "semi-abs", max_nodes = NULL)
    }
    if ("aver-abs-Wil-value" %in% methods_unpair) {
      rejection["aver-abs-Wil-value"] = 
        sudo_wilcoxon(data.frame(dat), n_permute = n_permute,
                      rank_type = FALSE, alg_type = "semi-abs", max_nodes = NULL)
    }}
    if (any(grepl("res", methods_unpair))) {
      df = data.frame(Y = e, Z = dat$Z, X = dat$X)
      p_res_vector = batch_pseudo_wilcoxon(dat = df, alg_type = alg_type, n_permute = n_permute)
      for(avail_method in names(p_res_vector)){
        if (paste(avail_method, "_res", sep = "") %in% methods_unpair) {
          rejection[paste(avail_method, "_res", sep = "")] =
            p_res_vector[avail_method]
        }
      }
    }
    if (any(grepl("outcome", methods_unpair))) {
      df = data.frame(Y = dat$Y, Z = dat$Z, X = dat$X)
      p_outcome_vector = batch_pseudo_wilcoxon(df, n_permute = n_permute)
      for(avail_method in names(p_outcome_vector)){
        if (paste(avail_method, "_outcome", sep = "") %in% methods_unpair) {
          rejection[paste(avail_method, "_outcome", sep = "")] =
            p_outcome_vector[avail_method]
        }
      }
    }
    
    return(rejection)
  }
  rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  # rejections = lapply(1:R, wrapper_func)
  return(rejections)
}



experiment_mixed = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    pvals = vector(length = length(methods_mixed)); names(pvals) = methods_mixed
    
    dat_pair = sample_pair_side(t = t, n = frac_pair*n, m = frac_pair*m, delta = delta,
                                prob = prob_assign, C_delta = C_delta, e_D = e_D, dat_type = dat_type)
    dat_merge = list(Z = c(dat_pair$Z1, dat_pair$Z2), Y = c(dat_pair$Y1, dat_pair$Y2),
                     X = rbind(dat_pair$X1, dat_pair$X2))
    if (res_type == "linear") {
      df = data.frame(Y = dat_pair$Y1 - dat_pair$Y2, X = dat_pair$X1 - dat_pair$X2)
      e_pair = lm(Y ~ .^2 -1, data = df)$residuals
    } else if (res_type == "nonlinear") {
      df = data.frame(Y = dat_pair$Y1, X = dat_pair$X1)
      e1 = df$Y- randomForest(Y ~., data = df)$predicted
      df = data.frame(Y = dat_pair$Y2, X = dat_pair$X2)
      e2 = df$Y - randomForest(Y ~., data = df)$predicted
      e_pair = e1 - e2
    }
    dat_unpair = sample_unpair(t = t, n = ceiling((1 - frac_pair)*n), m = (1 - frac_pair)*m, delta = delta,
                        prob = prob_assign, C_delta = C_delta, e_D = e_D, dat_type = dat_type)
    df = data.frame(Y = dat_unpair$Y, dat_unpair$X);
    if (res_type %in% c("linear")) {
      e_unpair = lm(Y ~ .^2 - 1, data = df)$residuals
    } else if(res_type == "nonlinear"){
      e_unpair = dat_unpair$Y - randomForest(Y ~ ., data = df)$predicted
    } 
    
    dat_sep = list(Z = c(dat_merge$Z, dat_unpair$Z), Y = c(dat_merge$Y, dat_unpair$Y),
               X = c(dat_merge$X, dat_unpair$X))
    dat_diff = list(Z = c((dat_pair$Z1 - dat_pair$Z2 + 1)/2, dat_unpair$Z),
                    Y = c(dat_pair$Y1 - dat_pair$Y2, dat_unpair$Y),
                    e = c(e_pair, e_unpair))
    if ("MannW" %in% methods_mixed) {
      pvals["MannW"] = wilcox.test(dat_diff$Y ~ dat_diff$Z, paired = FALSE)$p.value
    }
    if ("Wilcoxon" %in% methods_mixed) {
      pvals["Wilcoxon"] = 
        wilcox.test(dat_diff$Y*(2*dat_diff$Z - 1), rep(0, n), paired = TRUE)$p.value
    }
    if ("covadj-MannW" %in% methods_mixed) {
      pvals["covadj-MannW"] = 
        wilcox.test(dat_diff$e ~ dat_diff$Z, paired = FALSE)$p.value
    }
    if ("covadj-Wilcoxon" %in% methods_mixed) {
      pvals["covadj-Wilcoxon"] = 
        wilcox.test((2*dat_diff$Z - 1)*dat_diff$e, rep(0, n), paired = TRUE)$p.value 
    }
    if ("covadj-Wilcoxon-value" %in% methods_mixed) {
      pvals["covadj-Wilcoxon-value"] = 
        permute_wilcoxon(Y_diff = dat_diff$e, Z_diff = (2*dat_diff$Z - 1), n_permute = n_permute) 
    }
    if ("indep-test" %in% methods_mixed) {
      pvals["indep-test"] = dcov.test(dat_diff$e, dat_diff$Z, R = 100)$p.value 
    }
    if ("cor-test" %in% methods_mixed) {
      pvals["cor-test"] = cor.test(dat_diff$e, dat_diff$Z)$p.value
    }
    if ("martingale-mixed" %in% methods_mixed) {
      pvals["martingale-mixed"] =
        split_mixed(dat_pair = dat_merge, dat_unpair = dat_unpair, ub = ub)
    }
    if ("pseudo-Wilcoxon-rank" %in% methods_mixed) {
      pvals["pseudo-Wilcoxon-rank"] =
        sudo_wilcoxon_mixed(data.frame(dat_sep), n_permute = n_permute, n_pair = length(dat_pair$Z1),
                      rank_type = TRUE, alg_type = "sudo", max_nodes = NULL) 
    }
    if ("pseudo-Wilcoxon-value" %in% methods_mixed) {
      pvals["pseudo-Wilcoxon-value"] = 
        sudo_wilcoxon_mixed(data.frame(dat_sep), n_permute = n_permute, n_pair = length(dat_pair$Z1),
                      rank_type = FALSE, alg_type = "sudo", max_nodes = NULL)
    }
    if ("interactive" %in% methods_mixed) {
      pvals["interactive"] = 
        i_mixed(dat_pair = dat_pair, dat_unpair = dat_unpair,
                ub = ub, alternative = "two_sided", prob_assign = prob_assign)
    }
    return(pvals)
  }
  #rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  pvals = lapply(1:R, wrapper_func)
  return(pvals)
}


experiment_pval = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  rejection = matrix(ncol = 3, nrow = R); 
  colnames(rejection) = c("i-Mann-Whitney", "i-Mann-Whitney-pval", "order-Mann-Whitney-pval")
  for (i in 1:R) {
    print(i)
    dat = sample_unpair(t = t, n = n, m = m, delta = delta, C_delta = C_delta, C_f = C_f, e_D = e_D,
                        blind = blind)
    rejection[i, "i-Mann-Whitney"] = i_mann(X = dat$X, Z = dat$Z, Y = dat$Y, ub = ub)
    rejection[i, "i-Mann-Whitney-pval"] = i_mann_pval(X = dat$X, Z = dat$Z, Y = dat$Y) < alpha
    rejection[i, "order-Mann-Whitney-pval"] = order_mann_pval(Z = dat$Z, Y = dat$Y) < alpha
  }
  return(rejection)
}

experiment_multi = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  rejection = matrix(ncol = 2, nrow = R); 
  colnames(rejection) = c("Friedman", "i-Friedman")
  for (i in 1:R) {
    print(i)
    dat = sample_multi_pair(n = n, p = p, delta = delta, e_D = e_D)
    
    rejection[i, "Friedman"] = friedman.test(as.vector(dat$Y), groups = as.vector(dat$Z),
                                             blocks = rep(1:n, 3))$p.value < alpha
    
    rejection[i, "i-Friedman"] = i_friedman(Y_mat = dat$Y, Z_mat = dat$Z, ub = ub)
  }
  return(rejection)
}





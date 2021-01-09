# repeat experiments for interactive methods
# para_vary: a list of parameters that are different from the default in input.R file
experiment_interactive = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    p_val = vector(length = length(methods_interactive)); names(p_val) = methods_interactive
    
    dat = sample_generator(n = n, m = m,
                        treatment_type = treatment_type, control_type = control_type, 
                        C_delta = C_delta, e_D = e_D, 
                        redundant = redundant, blind = blind)
    
    if ("CovAdj-Wilcoxon-linear" %in% methods_interactive) {
      e = lm(Y ~ (. - A)^2 - 1, data = dat)$residuals
      p_val["CovAdj-Wilcoxon-linear"] = wilcox.test(e ~ dat$A, paired = FALSE)$p.value 
    }
    if ("CovAdj-Wilcoxon-robust" %in% methods_interactive) {
      e = rlm(Y ~ (. - A)^2 - 1, data = dat)$residuals
      p_val["CovAdj-Wilcoxon-robust"] = wilcox.test(e ~ dat$A, paired = FALSE)$p.value 
    }
    if ("CovAdj-Wilcoxon-quadratic" %in% methods_interactive) {
      e = rlm(Y ~ (. - A)^2 + I(X.3^2) - 1, data = dat)$residuals
      p_val["CovAdj-Wilcoxon-quadratic"] = wilcox.test(e ~ dat$A, paired = FALSE)$p.value 
    }
    
    if ("linear-CATE-test" %in% methods_interactive) {
      p_val["linear-CATE-test"] = linear_CATE_test(dat)
    }
    
    if ("i-Wilcoxon-linear" %in% methods_interactive) {
      p_val["i-Wilcoxon-linear"] = 
        i_Wilcoxon(dat, alg_type = "linear", iter_round = iter_round)
    }
    if ("i-Wilcoxon-robust" %in% methods_interactive) {
      p_val["i-Wilcoxon-robust"] = 
        i_Wilcoxon(dat, alg_type = "robust", iter_round = iter_round)
    }
    if ("i-Wilcoxon-quadratic" %in% methods_interactive) {
      p_val["i-Wilcoxon-quadratic"] = 
        i_Wilcoxon(dat, alg_type = "quadratic", iter_round = iter_round)
    }
    
    if ("i-Wilcoxon-linear-signedA" %in% methods_interactive) {
      p_val["i-Wilcoxon-linear-signedA"] = 
        i_Wilcoxon(dat, alg_type = "linear", sum_type = "signed_A", order_by = "abs_pred",
                   iter_round = iter_round)
    }
    if ("i-Wilcoxon-robust-signedA" %in% methods_interactive) {
      p_val["i-Wilcoxon-robust-signedA"] = 
        i_Wilcoxon(dat, alg_type = "robust", sum_type = "signed_A", order_by = "abs_pred", 
                   iter_round = iter_round)
    }
    if ("i-Wilcoxon-robust-signedA_adapt" %in% methods_interactive) {
      p_val["i-Wilcoxon-robust-signedA_adapt"] = 
        i_Wilcoxon(dat, alg_type = "robust", sum_type = "signed_A", order_by = "adapt_pred",
                   iter_round = iter_round)
    }
    if ("i-Wilcoxon-quadratic-signedA" %in% methods_interactive) {
      p_val["i-Wilcoxon-quadratic-signedA"] = 
        i_Wilcoxon(dat, alg_type = "quadratic", sum_type = "signed_A", order_by = "abs_pred",
                   iter_round = iter_round)
    }
    if ("i-Wilcoxon-quadratic-signedA_adapt" %in% methods_interactive) {
      p_val["i-Wilcoxon-quadratic-signedA_adapt"] = 
        i_Wilcoxon(dat, alg_type = "quadratic", sum_type = "signed_A", order_by = "adapt_pred",
                   iter_round = iter_round)
    }
    
    if ("i-Wilcoxon-oracle" %in% methods_interactive) {
      p_val["i-Wilcoxon-oracle"] = 
        i_Wilcoxon(dat, C_delta = C_delta, alg_type = "oracle",
                   sum_type = "signed_A", order_by = "abs_pred")
    }
    return(p_val)
  }
  p_vals = lapply(1:R, wrapper_func)
  # p_vals = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  return(p_vals)
}


# repeat experiments for non-interactive variants of the Wilcoxon signed rank tests
# para_vary: a list of parameters that are different from the default in input.R file
experiment_var_Wilcoxon = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    dat = sample_generator(n = n, m = m,
                        treatment_type = treatment_type, control_type = control_type, 
                        C_delta = C_delta, e_D = e_D, 
                        redundant = redundant, blind = blind)
    if (alg_type == "linear") {
      R = rlm(Y ~ (. - A)^2 - 1, data = dat)$residuals
    } else if(alg_type == "RF"){
      R = dat$Y - randomForest(Y ~ . - A, data = dat)$predicted
    } 
    
    df = cbind(dat, R = R)
    p_vals = var_wilcoxon(dat = df, alg_type = alg_type, n_permute = n_permute)
    #return(p_vals[methods_var_Wilcoxon])
    return(p_vals)
  }
  p_vals = lapply(1:R, wrapper_func)
  return(p_vals)
}






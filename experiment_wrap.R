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
        i_Wilcoxon(dat, alg_type = "linear")
    }
    if ("i-Wilcoxon-robust" %in% methods_interactive) {
      p_val["i-Wilcoxon-robust"] = 
        i_Wilcoxon(dat, alg_type = "robust")
    }
    if ("i-Wilcoxon-quadratic" %in% methods_interactive) {
      p_val["i-Wilcoxon-quadratic"] = 
        i_Wilcoxon(dat, alg_type = "quadratic")
    }
    return(p_val)
  }
  p_vals = mclapply(1:R, wrapper_func, mc.cores = detectCores())
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
      R = lm(Y ~ (. - A)^2 - 1, data = dat)$residuals
    } else if(alg_type == "RF"){
      R = dat$Y - randomForest(Y ~ . - A, data = dat)$predicted
    } 
    
    df = cbind(dat, R = R)
    p_vals = var_wilcoxon(dat = df, alg_type = alg_type, n_permute = n_permute)
    return(p_vals[methods_var_Wilcoxon])
  }
  p_vals = lapply(1:R, wrapper_func)
  return(p_vals)
}






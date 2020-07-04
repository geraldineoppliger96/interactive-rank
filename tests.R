##############################################################
#################### interactive Wilcoxon ####################
##############################################################
# i-Wilcoxon
# dat:         a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# alg_type:    type of algorithm for modeling and ordering
# alterantive: one-sided or two-sided boundary (default as two-sided)
i_Wilcoxon = function(dat, alg_type = "robust", alternative = "two_sided"){
  n = length(dat$A)
  if(alg_type == "linear") {
    pred_func = em_linear
  } else if (alg_type == "robust"){
    pred_func = em_robust
  } else if (alg_type == "quadratic"){
    pred_func = em_quadratic
  }
  
  cand_set <- rep(FALSE, n); ordered_A <- vector(length = n); k <- 1
  while(k <= 500){
    if (k %% 100 == 1) {
      pred <- pred_func(dat, cand_set = cand_set, iter = 20)
      #pred = pred - min(pred) + 1 #make sure pred score strictly > 0!!!!!!
    }
    inc_ind <- which.max(pred*(!cand_set) + (-1)*cand_set); cand_set[inc_ind] <- TRUE
    ordered_A[k] <- dat$A[inc_ind]
    k = k + 1
  }
  k = 1:n; S_cum = cumsum(2*ordered_A - 1)
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(2/n)*k + sqrt(n/8)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(2/n)*k + sqrt(n/8)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(2/n)*k + sqrt(n/8)))^2))
  }
  return(p_val)
}


# em_linear: EM algorithm for a regular linear model
# dat:      a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set: a vector of indicators for whether subject i is ordered 
# iter:     number of interations within a single EM algorithm
em_linear <- function(dat, cand_set, iter){
  n <- length(dat$A); Y <- dat$Y - min(dat$Y)
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * dnorm(dat$Y - mu_1)/(dnorm(dat$Y - mu_1) + dnorm(dat$Y - mu_0)) +
      cand_set*dat$A
    w[is.na(w)] = 1/2
    if(any(w != 1)){
      mu_0 = lm(Y ~ (. - A)^2, data = dat, weights = 1 - w)$fitted.values
    }
    mu_1 = lm(Y ~ (. - A)^2, data = dat, weights = w)$fitted.values
  }
  return(w)
}


# em_robust: EM algorithm for a robust linear model
# dat:      a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set: a vector of indicators for whether subject i is ordered 
# iter:     number of interations within a single EM algorithm
em_robust <- function(dat, cand_set, iter){
  n <- length(dat$A); Y <- dat$Y - min(dat$Y); 
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * dnorm(dat$Y - mu_1)/(dnorm(dat$Y - mu_1) + dnorm(dat$Y - mu_0)) +
      cand_set*dat$A
    w[is.na(w)] = 1/2
    if (any(w != 1)) {
      mu_0 = rlm(Y ~ (. - A)^2, data = dat, weights = 1 - w)$fitted.values
    }
    mu_1 = rlm(Y ~ (. - A)^2, data = dat, weights = w)$fitted.values
  }
  return(w)
}


# em_quadratic: EM algorithm for a robust model with quadratic term for X(3)
# dat:      a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set: a vector of indicators for whether subject i is ordered 
# iter:     number of interations within a single EM algorithm
em_quadratic <- function(dat, cand_set, iter){
  n <- length(dat$A); Y <- dat$Y - min(dat$Y); 
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * dnorm(dat$Y - mu_1)/(dnorm(dat$Y - mu_1) + dnorm(dat$Y - mu_0)) +
      cand_set*dat$A
    w[is.na(w)] = 1/2
    if(any (w != 1)){
      mu_0 = rlm(Y ~ (. - A)^2 + I(X.3^2), data = dat, weights = 1 - w)$fitted.values
      mu_0[is.na(mu_0)] = 0
    }
    mu_1 = rlm(Y ~ (. - A)^2 + I(X.3^2), data = dat, weights = w)$fitted.values
  }
  return(w)
}



##############################################################
#################### linear-CATE-test ########################
##############################################################
linear_CATE_test = function(dat){
  # dat: a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
  n = length(dat$A)
  
  Y_hat <- lm(Y ~ (. - A)^2 , data = dat)$fitted.values 
  X_cov = model.matrix(Y ~ (. - A)^2, data = dat)
  
  B <- diag(as.vector((dat$A - 1/2)*(dat$Y - Y_hat))) %*% X_cov; b <- colMeans(B)
  chi_stat <- tryCatch({b %*% solve(cov(B)/n) %*% b}, error = function(cond){return(NA)})
  p_val = pchisq(chi_stat, df = ncol(X_cov), lower.tail = FALSE)
  return(p_val)
}



##############################################################
#################### variants of Wilcoxon ####################
##############################################################
# var-Wilcoxon: non-interactive variants of the Wilcoxon signed-rank test
# dat:       a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# alg_type:  type of algorithm for modeling (default as random forest)
# n_permute: number of permutations to get a single p-value (default as 200)
var_wilcoxon = function(dat, alg_type = "RF", n_permute = 200) {
  W_obs = stat_generator(dat, alg_type = alg_type); print(W_obs)
  W_permute = foreach(i = 1:n_permute, .combine = cbind,
                      .export = c("stat_generator"),
                      .packages = c("foreach", "randomForest")) %dopar% {
                        dat_permute = dat; 
                        dat_permute$A = base::sample(dat_permute$A)
                        stat_generator(dat_permute, alg_type = alg_type)
                      }
  pval = rowMeans(apply(W_permute, 2, function(x) {x > W_obs}))
  return(pval)
}


# stat_generator: compute test statstistic W using several choices of E_i,
#                 given the observed or permuted dataset
# dat:       a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# alg_type:  type of algorithm for modeling (default as random forest)
stat_generator = function(dat, alg_type) {
  W = vector(length = 5)
  names(W) = c("R(X)",
               "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
               "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
               "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
  
  if (alg_type == "linear") {
    model_res = lm(R ~ . - Y, data = dat)
    model_outcome = lm(Y ~ . - R, data = dat)
  } else if (alg_type == "RF") {
    model_res = randomForest(R ~ . - Y, data = dat)
    model_outcome = randomForest(Y ~ . - R, data = dat)
  }
  rhat_true = predict(model_res); yhat_true = predict(model_outcome)
  dat_false = dat; dat_false$A = 1 - dat_false$A
  rhat_false = predict(model_res, newdata = dat_false)
  yhat_false = predict(model_outcome, newdata = dat_false)
  
  E = (2*dat$A - 1)*dat$R; W["R(X)"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  
  E = (2*dat$A - 1)*(dat$Y - yhat_false)
  W["R(X, 1-A)"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  
  E = (2*dat$A - 1)*(dat$R - rhat_false)
  W["R - hat(R)(X, 1 - A)"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  
  E = abs(rhat_false - dat$R) - abs(rhat_true - dat$R)
  W["|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  
  comb_sign = (2*dat$A - 1)*(dat$R - rhat_false) >= 0 |
    (2*dat$A - 1)*(rhat_true - rhat_false) >= 0
  E = (2*comb_sign - 1)*(abs(rhat_false - dat$R) - abs(rhat_true - dat$R))
  W["S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  
  return(W)
}



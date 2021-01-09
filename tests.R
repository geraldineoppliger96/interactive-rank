##############################################################
#################### interactive Wilcoxon ####################
##############################################################
# i-Wilcoxon
# dat:         a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# alg_type:    type of algorithm for modeling and ordering
# order_by:    strategy for ordering
# sum_type:    strategy to decide weights
# C_delta:     scale of treatment effect, used when alg_type = "oracle"
# iter_round:  number of iterations we skip to update the modeling
i_Wilcoxon = function(dat, alg_type = "robust", order_by = "pred", sum_type = "A",
                      C_delta = NA, iter_round = 100){
  n = length(dat$A)
  if(alg_type == "linear") {
    pred_func = em_linear
  } else if (alg_type == "robust") {
    pred_func = em_robust
  } else if (alg_type == "quadratic") {
    pred_func = em_quadratic
  } else if (alg_type == "oracle") {
    pred_func = em_oracle
  }
  
  cand_set <- rep(FALSE, n); ordered_A <- vector(length = n); ordered_pred <- vector(length = n)
  k <- 1; X_seq = vector(length = 0); S_cum = rep(0, n); increase = TRUE; model_correct = TRUE
  while(k <= n){
    if (k %% min(iter_round, n/5) == 1) {
      pred <- pred_func(dat, cand_set = cand_set, C_delta = C_delta, iter = 20)
    }
    if (order_by == "pred") {
      inc_ind <- which.max(pred*(!cand_set) + (-1)*cand_set); cand_set[inc_ind] <- TRUE
      ordered_A[k] <- dat$A[inc_ind]
      ordered_pred[k] <- pred[inc_ind]
    } else if (order_by == "abs_pred") {
      inc_ind <- which.max(abs(pred - 1/2)*(!cand_set) + (-1)*cand_set); cand_set[inc_ind] <- TRUE
      ordered_A[k] <- dat$A[inc_ind]
      ordered_pred[k] <- pred[inc_ind]
    } else if (order_by == "adapt_pred") {
      inc_ind <- which.max(abs(pred - 1/2)*(!cand_set) + (-1)*cand_set); cand_set[inc_ind] <- TRUE
      if(k > iter_round & mean(X_seq[1:(k-1)]) < 0) {
        ordered_A[k] <- dat$A[inc_ind]
        ordered_pred[k] <- 1 - pred[inc_ind]
      } else {
        ordered_A[k] <- dat$A[inc_ind]
        ordered_pred[k] <- pred[inc_ind]
      }
    }
    X_seq = c(X_seq, (2*ordered_A[k] - 1)*(2*(ordered_pred[k] > 1/2) - 1))
    k = k + 1
  }
  k = 1:n
  if (sum_type == "A") {
    X_seq = 2*ordered_A - 1
  } else if (sum_type == "signed_A") {
    X_seq = (2*ordered_A - 1)*(2*(ordered_pred > 1/2) - 1)
  }
  
  S_cum = cumsum(X_seq)
  p_val = 2*min(exp(-(abs(S_cum) / (sqrt(2/n)*k + sqrt(n/8)))^2))
  return(p_val)
}


# em_linear: EM algorithm for a regular linear model
# dat:      a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set: a vector of indicators for whether subject i is ordered 
# iter:     number of interations within a single EM algorithm
em_linear <- function(dat, cand_set, iter, C_delta = NA){
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
em_robust <- function(dat, cand_set, iter, C_delta = NA){
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
  w[is.na(w)] = 1/2
  return(w)
}

# em_oracle: guess the treatment assignment with oracle knowledge of potential outcomes
# dat:      a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set: a vector of indicators for whether subject i is ordered 
# iter:     number of interations within a single EM algorithm
# C_delta:  scale of the treatment effect
em_oracle <- function(dat, cand_set, iter, C_delta){
  mu_0 = 5*(dat$X.1 + dat$X.2 + dat$X.3)
  mu_1 = C_delta*(dat$X.1*dat$X.2 + dat$X.3)*2/5 + mu_0
  w = (!cand_set) * dnorm(dat$Y - mu_1)/(dnorm(dat$Y - mu_1) + dnorm(dat$Y - mu_0)) +
    cand_set*dat$A
  w[is.na(w)] = 1*(dat$Y[is.na(w)] > 0)
  return(w)
}

# em_quadratic: EM algorithm for a robust model with quadratic term for X(3)
# dat:      a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set: a vector of indicators for whether subject i is ordered 
# iter:     number of interations within a single EM algorithm
em_quadratic <- function(dat, cand_set, iter, C_delta = NA){
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
  
  Y_hat <- lm(Y ~ (. - A)^2, data = dat)$fitted.values 
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
                      .packages = c("foreach", "randomForest", "MASS")) %dopar% {
                        dat_permute = dat; 
                        dat_permute$A = base::sample(dat_permute$A)
                        stat_generator(dat_permute, alg_type = alg_type)
                      }
  pval = rowMeans(apply(W_permute, 2, function(x) {x > W_obs}), na.rm = TRUE)
  return(pval)
}


# stat_generator: compute test statstistic W using several choices of E_i,
#                 given the observed or permuted dataset
# dat:       a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# alg_type:  type of algorithm for modeling (default as random forest)
stat_generator = function(dat, alg_type) {
  W = rep(NA, 5)
  names(W) = c("R(X)",
               "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
               "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
               "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
  
  if (alg_type == "linear") {
    model_res = tryCatch({rlm(R ~ (. - Y)^2, data = dat)}, error = function(cond){return(NA)})
    model_outcome = tryCatch({rlm(Y ~ (. - R)^2, data = dat)}, error = function(cond){return(NA)})
  } else if (alg_type == "RF") {
    model_res = randomForest(R ~ . - Y, data = dat)
    model_outcome = randomForest(Y ~ . - R, data = dat)
  }
  if(!is.na(model_res) & !is.na(model_outcome)) {
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
  } 
  return(W)
}



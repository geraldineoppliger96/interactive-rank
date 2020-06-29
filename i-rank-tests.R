i_sign_rank = function(Y1, Y2, ub, alpha, x = NA){
  S = 2*(Y1 > Y2) - 1
  M = abs(Y1 - Y2)
  
  abs_sum = S[order(M, decreasing = TRUE)] %>% cumsum %>% abs
  reject = any(abs_sum > ub)
  return(reject)
}

i_wilcoxon = function(X1, Z1, Y1, X2, Z2, Y2, ub, alpha, dat_type, alg_type){
  n = length(Z1)
  
  df = data.frame(Y = Y1 - Y2, X = X1 - X2)
  if (dat_type == "linear") {
    e = lm(Y~ .^2 - 1, data = df)$residuals
  } else if (dat_type == "nonlinear") {
    e = df$Y - randomForest(Y ~., data = df)$predicted
  }
  
  if (alg_type == "em-lm") {
    pred_func = em_diff
  } else if (alg_type == "em-rf") {
    pred_func = em_rf
  }
  
  d = (Z1 - Z2)*e; S = 2*(d > 0) - 1

  rej = FALSE; cand_set = rep(FALSE, n); sign_sum = 0; k = 1
  while(!rej & k <= n){
    d_obs <- abs(d)*(!cand_set) + d*(cand_set)
    if (k %% 100 == 1) {
      pred <- pred_func(d_obs, (X1 + X2)/2, cand_set, iter = 10)
      pred = pred - min(pred) + 1
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    sign_sum = sign_sum + S[inc_ind]
    rej <- (abs(sign_sum) > ub[k])
    k = k + 1
  }
  return(rej)
}

em_diff <- function(d_obs, X, cand_set, iter){
  mu = rep(0.1, length(d_obs))
  for (i in 1:iter) {
    w = 1/(1 + exp(-2*mu*d_obs))
    w[cand_set] <- 1 # it should be the sign!!!!!!!
    #mu = lm((2*w - 1)*d_obs ~ (X[,1] + X[,2] + X[,3])^2)$fitted.values
    weight_df = data.frame(Y = (2*w - 1)*d_obs, X = X)
    mu = lm(Y ~ .^2, data = weight_df)$fitted.values
  }
  return(w)
}

em_rf <- function(d_obs, X, cand_set, iter){
  n = length(d_obs)
  prob = runif(n)
  
  predictor = cbind(d_obs, X)
  reg = randomForest(y = prob, x = predictor)   
  prob = predict(reg)
  prob[is.na(prob)] = 1/2
}


## predict assignment
i_wilcoxon_two = function(X1, Z1, Y1, X2, Z2, Y2, ub, alpha, prob_assign){
  S = (Z1 - Z2 + 1)/2; S = (S - prob_assign)/sqrt(prob_assign*(1 - prob_assign))
  n = length(S); cand_set = rep(FALSE, n)
  rej = FALSE; sign_sum = 0; k = 1
  while(!rej & k <= n){
    if (k %% 100 == 1) {
      pred <- em_pair_sep(X1, Z1, Y1, X2, Z2, Y2, cand_set, iter = 10)
      pred = pred - min(pred) + 1
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    sign_sum = sign_sum + S[inc_ind]
    rej <- (abs(sign_sum) > ub[k])
    k = k + 1
  }
  return(rej)
}

em_pair_both <- function(X1, Z1, Y1, X2, Z2, Y2, cand_set, iter){
  n <- length(Y1)
  mu_1 = rep(0.1, 2*n); mu_0 = rep(0, 2*n)
  for (i in 1:iter) {
    f_p = exp((Y1 - mu_1[1:n])^2/2) * exp((Y2 - mu_0[(n + 1):(2*n)])^2/2)
    f_n = exp((Y1 - mu_0[1:n])^2/2) * exp((Y2 - mu_1[(n + 1):(2*n)])^2/2)
    w = (!cand_set) * (f_p / (f_p + f_n)) + cand_set*(Z1 - Z2)
    
    df = data.frame(Y = c(Y1, Y2), X = rbind(X1, X2))
    mu_1 = lm(Y ~ .^2, data = df,
              weights = c(w, 1 - w))$fitted.values
    mu_0 = lm(Y ~ .^2, data = df,
              weights = c(1 - w, w))$fitted.values
    
  }
  #print(unique(mu))
  return(w)
}

em_pair_sep <- function(X1, Z1, Y1, X2, Z2, Y2, cand_set, iter){
  Y = c(Y1, Y2); X = rbind(X1, X2); Z = c(Z1, Z2); cand_set = rep(cand_set, 2)
  
  n <- length(Y); Y <- Y - min(Y)
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * 1/(1 + exp((mu_0 - mu_1)*Y + mu_1^2/2 - mu_0^2/2)) + cand_set*Z
    df = data.frame(Y, X)
    mu_1 = lm(Y ~ .^2, data = df,
              weights = w)$fitted.values
    if(any(w != 1)){
      mu_0 = lm(Y ~ .^2, data = df,
                weights = 1 - w)$fitted.values
    }
  }
  
  w_pair = w[1:(n/2)]/(w[1:(n/2)] + w[(n/2 + 1):n])
  w_pair[is.na(w_pair)] = 0.5
  #print(unique(mu))
  return(w_pair)
}


## split
i_wilcoxon_split = function(dat, ub){
  D = as.data.frame(dat)
  n = nrow(D); split_ind = sample.split(1:n, SplitRatio =  0.5)
  D = add_column(D, Z = as.factor(D$Z1 - D$Z2)); D = subset(D, select = -c(Z1, Z2))
  D1 = subset(D, split_ind == TRUE); D2 = subset(D, split_ind == FALSE)
  
  n2 = nrow(D2); S = 2*as.numeric(D2$Z) - 3
  rej = FALSE; cand_set = rep(FALSE, n2); sign_sum = 0; k = 1
  while(!rej & k <= n2){
    if (k %% 100 == 1) {
      D1 = rbind(D1, subset(D2, cand_set == TRUE))
      classifier = randomForest(Z ~., data = D1)
      pred = predict(classifier, newdata = D2, type = "vote")[,2]
      pred = pred - min(pred) + 1
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    sign_sum = sign_sum + S[inc_ind]
    rej <- (abs(sign_sum) > ub[k])
    k = k + 1
  }
  return(rej)
}



######### unpaired #################
i_mann = function(X, Z, Y, ub, type, oracle){
  n = length(Z)
  if(type == "cont") {
    pred_func = em_unpair_both
  } else if (type == "bi"){
    pred_func = em_unpair_bi_both
  }
  
  rej = FALSE; cand_set <- rep(FALSE, n); sign_sum <- 0; k <- 1
  while(!rej & k <= n){
    if (k %% 100 == 1) {
      pred <- pred_func(Y = Y, Z = Z, X = X, cand_set = cand_set, iter = 5, oracle = oracle)
      pred = pred - min(pred) + 1
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    sign_sum <- sign_sum + (2*Z[inc_ind] - 1)
    rej = abs(sign_sum) > ub[k]
    k = k + 1
  }
  return(rej)
}


em_unpair_both <- function(Y, Z, X, cand_set, iter, oracle){
  n <- length(Y); Y <- Y - min(Y); 
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * 1/(1 + exp((mu_0 - mu_1)*Y + mu_1^2/2 - mu_0^2/2)) + cand_set*Z
    w[is.na(w)] = 1/2
    if(any (w != 1)){
      df_0 = data.frame(Y, X)
      mu_0 = rlm(Y ~ .^2 + I(X3^2), data = df_0,
                weights = 1 - ((!cand_set) * w + cand_set*Z))$fitted.values
      mu_0[is.na(mu_0)] = 0
      # mu_0 = glm(Y ~ .^2, data = df, family = poisson,
      #           weights = 1 - ((!cand_set) * w + cand_set*Z))$fitted.values
    }
    # res = Y - mu_0;
    # X_1 = X; X_1[,3] = log(X_1[,3] - min(X_1[,3]) + 1)
    # df_1 = data.frame(Y = log(res - min(res) + 1), X_1)
    # mu_1 = rlm(Y ~ .^2, data = df_1,
    #            weights = (!cand_set) * w + cand_set*Z)$fitted.values
    # mu_1 = mu_0 + exp(mu_1) + min(res) - 1
    df_1 = data.frame(Y = Y, X)
    mu_1 = rlm(Y ~ .^2 + I(X3^2), data = df_1, weights = (!cand_set) * w + cand_set*Z)$fitted.values
    # df_1 = data.frame(Y = res, X)
    # mu_1 = lm(Y ~ .^2, data = df_1,
    #           weights = (!cand_set) * w + cand_set*Z)$fitted.values
    # mu_1 = mu_0 + mu_1
  }
  if (0) {
    par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
    p = plot(lm(Y ~ .^2, data = df_0,
            weights = 1 - ((!cand_set) * w + cand_set*Z)), 5, caption = "")
    ggsave(filename = paste(dirname(getwd()),"/figure/interactive_qq.png", sep = ""),
           plot = p, width = 4, height = 3.6)

    plot(lm(Y ~ .^2, data = df_1,
             weights = (!cand_set) * w + cand_set*Z), 5, caption = "")
    par(mfrow=c(1,1))
    
    plot(mu_0, Y - mu_0, xlab = "Fitted values", ylab = "Residuals",
         ylim = c(-5, 5), xlim = c(min(mu_1), max(mu_1)))
    df_res = data.frame(y = Y - mu_1, x = mu_1)
    loess_fit <- loess(y ~ x, df_res)
    lines(sort(df_res$x), predict(loess_fit)[order(df_res$x)], col = "red", lwd = 2)
  }
  
  
  # d = ncol(X); n_para = d*(d-1) + 2*d + 1
  # BIC = -2*sum(log(w*dnorm(Y - mu_1) + (1 - w)*dnorm(Y - mu_0))) + n_para*log(n)
  #print(unique(mu))
  return(w)
}


em_unpair_both_default <- function(Y, Z, X, cand_set, iter, oracle){
  n <- length(Y); Y <- Y - min(Y); 
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * 1/(1 + exp((mu_0 - mu_1)*Y + mu_1^2/2 - mu_0^2/2)) + cand_set*Z
    w[is.na(w)] = 1/2
    if(any(w != 1)){
      df_0 = data.frame(Y, X)
      mu_0 = rlm(Y ~ .^2, data = df_0,
                weights = 1 - ((!cand_set) * w + cand_set*Z))$fitted.values
    }
    df_1 = data.frame(Y = Y, X)
    mu_1 = rlm(Y ~ .^2, data = df_1,
              weights = (!cand_set) * w + cand_set*Z)$fitted.values
  }
  return(w)
}

em_unpair_both_naive <- function(Y, Z, X, cand_set, iter, oracle){
  n <- length(Y); Y <- Y - min(Y); 
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  for (i in 1:iter) {
    w = (!cand_set) * 1/(1 + exp((mu_0 - mu_1)*Y + mu_1^2/2 - mu_0^2/2)) + cand_set*Z
    w[is.na(w)] = 1/2
    if(any(w != 1)){
      df_0 = data.frame(Y, X)
      mu_0 = lm(Y ~ .^2, data = df_0,
                 weights = 1 - ((!cand_set) * w + cand_set*Z))$fitted.values
    }
    df_1 = data.frame(Y = Y, X)
    mu_1 = lm(Y ~ .^2, data = df_1,
               weights = (!cand_set) * w + cand_set*Z)$fitted.values
  }
  return(w)
}


em_unpair_bi_both <- function(Y, Z, X, cand_set, iter, oracle){ #y = 0/1
  n <- length(Y)
  if (oracle) {
    if(ncol(X) == 2){ X = cbind(X, rep(0, n))}
    mu_1 = 1/(1 + exp(-out_fun(t, X, 1)))
    mu_0 = 1/(1 + exp(-out_fun(t, X, 0)))
    w = (!cand_set) * 1/(1 + ((mu_0/mu_1)^Y * ((1 - mu_0)/(1 - mu_1))^(1 - Y))) + cand_set*Z
  } else {
    mu_1 = rep(0.6, n); mu_0 = rep(0.5, n)
    for (i in 1:iter) {
      w = (!cand_set) * 1/(1 + ((mu_0/mu_1)^Y * ((1 - mu_0)/(1 - mu_1))^(1 - Y))) + cand_set*Z

      df = data.frame(Y, X)
      mu_1 = glm(Y ~ .^2, data = df, family = binomial(),
                weights = w)$fitted.values
      if(any (w != 1)){
        mu_0 = glm(Y ~ .^2, data = df, family = binomial(),
                  weights = 1 - w)$fitted.values
      }
    }
  }
  #mu_1 = 1/(1 + exp(-C_delta*(X[,1]*X[,2] + 0.5*X[,3])))
  #mu_0 = 1/2

  return(w)
}


em_unpair_nonpara <- function(Y, Z, X, cand_set, iter, oracle){
  n <- length(Y); Y <- Y - min(Y)
  mu_1 = rep(0.1, n); mu_0 = rep(0, n)
  
  w_1 = (!cand_set) * 1/(1 + exp((mu_0 - mu_1)*Y + mu_1^2/2 - mu_0^2/2)) + cand_set*Z; w_0 = 1 - w_1
  df = data.frame(Y, X)
  model_1 = lm(Y ~ .^2, data = df, weights = w_1)
  if(any (w_1 != 1)){
    model_0 = lm(Y ~ .^2, data = df, weights = w_0)
    mu_0 = model_0$fitted.values
  }
  mu_1 = model_1$fitted.values 
  var_1 = summary(model_1)$sigma^2; var_0 = summary(model_0)$sigma^2
  ker_val = 1/h*outer(1:n, 1:n, Vectorize(function(i, j) dmvnorm((X[i, ] - X[j, ])/h)))  
  for (i in 1:iter) {
    r_1 = w_1*dnorm((Y - mu_1)/sqrt(var_1)) /
      (w_1*dnorm((Y - mu_1)/sqrt(var_1)) + w_0*dnorm((Y - mu_0)/sqrt(var_0)))
    r_0 = 1 - r_1
    
    a = sum(r_1 %*% ker_val); b = sum(r_0 %*% ker_val)
    w_1 = a/(a + b); w_0 = 1 - w_1
    mu_1 = sum((r_1 %*% ker_val) *Y)/a; mu_0 = sum((r_0 %*% ker_val)*Y)/b
    var_1 = sum((r_1 %*% ker_val)*(Y - mu_1))/a; var_0 = sum((r_0 %*% ker_val)*(Y - mu_0))/b
  }
  
  # d = ncol(X); n_para = d*(d-1) + 2*d + 1
  # BIC = -2*sum(log(w*dnorm(Y - mu_1) + (1 - w)*dnorm(Y - mu_0))) + n_para*log(n)
  #print(unique(mu))
  return(w_1)
}


CATE_test = function(Y, Z, X, alpha){
  n = length(Z)

  df = data.frame(Y, X)
  Y_hat <- lm(Y ~ .^2, data = df)$fitted.values 
  X_cov = model.matrix(Y ~ .^2, data = df)
  
  A <- diag(as.vector((Z - 1/2)*(Y - Y_hat))) %*% X_cov; a <- colMeans(A)
  chi_stat <- tryCatch({a %*% solve(cov(A)/n) %*% a}, error = function(cond){print(cond); return(NA)})
  rej = chi_stat > qchisq(1 - alpha, df = ncol(X_cov))
  
  return(rej)
}


Bonf_mann = function(Y, Z, X, ub, alpha, type = type){
  df = data.frame(Y, X)
  if(type == "cont"){
    e = lm(Y ~ .^2 - 1, data = df)$residuals
  } else {
    e = glm(Y ~ .^2 - 1, data = df, family = binomial())$residuals
  }
  
  rej_cov = wilcox.test(e ~ Z, paired = FALSE)$p.value < alpha/2
   
  rej_i = i_mann(X = X, Z = Z, Y = Y, ub = ub, type = type, oracle = FALSE)
  
  rej = (rej_cov | rej_i)
  return(rej)
}
  


####### multi-sample ########

i_friedman = function(Y_mat, Z_mat, ub){
  Y_1 = rowSums(Y_mat * (Z_mat == 1))
  Y_3 = rowSums(Y_mat * (Z_mat == 3))
  S = 2*(Y_1 > Y_3) - 1
  
  M = apply(Y_mat, 1, max) - apply(Y_mat, 1, min)
  
  abs_sum = S[order(M, decreasing = TRUE)] %>% cumsum %>% abs
  reject = any(abs_sum > ub)
  return(reject)
}



############# p-value; unpaired ##################

order_mann_pval = function(Z, Y, alternative = "two_sided"){
  n = length(Z); k = 1:n
  S_cum = cumsum(2*Z[order(Y, decreasing = TRUE)] - 1)
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(2/n)*k + sqrt(n/8)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(2/n)*k + sqrt(n/8)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(2/n)*k + sqrt(n/8)))^2))
  }
  return(p_val)
}

i_mann_pval = function(X, Z, Y, prob_assign, cut_prop = 1/4, alternative = "two_sided", type = "cont"){
  n = length(Z)
  if(type == "cont") {
    pred_func = em_unpair_both
  } else if (type == "bi"){
    pred_func = em_unpair_bi_both
  } else if (type == "default"){
    pred_func = em_unpair_both_default
  } else if (type == "naive"){
    pred_func = em_unpair_both_naive
  }

  cand_set <- rep(FALSE, n); ordered_Z <- vector(length = n); ordered_ind = vector(length = n); k <- 1
  while(k <= 500){
    if (k %% 100 == 1) {
      pred <- pred_func(Y = Y, Z = Z, X = X, cand_set = cand_set, iter = 20, oracle = FALSE)
      pred = pred - min(pred) + 1 #make sure pred score strictly > 0!!!!!!
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    ordered_Z[k] <- Z[inc_ind]
    ordered_ind[k] <- inc_ind
    k = k + 1
  }
  k = 1:n
  #S_cum = cumsum(2*ordered_Z - 1)
  S_cum = cumsum((ordered_Z - prob_assign)/sqrt(prob_assign*(1 - prob_assign)))
  m = n*cut_prop
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  }
  return(p_val)
}

i_mann_split_pval = function(X, Z, Y, cut_prop = 1/4, alternative = "two_sided"){
  D = data.frame(Z, Y, X); D[,1] = factor(D[,1])
  n = nrow(D); split_ind = sample.split(1:n, SplitRatio =  0.5)
  D1 = subset(D, split_ind == TRUE); D2 = subset(D, split_ind == FALSE)
  
  n = nrow(D2)
  cand_set <- rep(FALSE, n); ordered_Z <- vector(length = n); k <- 1
  while(k <= n){
    if (k %% 100 == 1) {
      D1 = rbind(D1, subset(D2, cand_set == TRUE))
      classifier = randomForest(Z ~., data = D1)
      pred = predict(classifier, newdata = D2, type = "vote")[,2]
      pred = pred - min(pred) + 1
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    ordered_Z[k] <- Z[inc_ind] 
    k = k + 1
  }
  k = 1:n
  S_cum = cumsum(2*ordered_Z - 1)
  m = n*cut_prop
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  }
  return(p_val)
}



CATE_test_pval = function(Y, Z, X){
  n = length(Z)
  
  df = data.frame(Y, X)
  Y_hat <- lm(Y ~ .^2, data = df)$fitted.values 
  X_cov = model.matrix(Y ~ .^2, data = df)
  
  A <- diag(as.vector((Z - 1/2)*(Y - Y_hat))) %*% X_cov; a <- colMeans(A)
  chi_stat <- tryCatch({a %*% solve(cov(A)/n) %*% a}, error = function(cond){return(NA)})
  #rej = chi_stat > qchisq(1 - alpha, df = ncol(X_cov))
  p_val = pchisq(chi_stat, df = ncol(X_cov), lower.tail = FALSE)
  return(p_val)
}


i_mixed = function(dat_pair, dat_unpair, ub, alternative = "two_sided", cut_prop = 1/4,
                   prob_assign = prob_assign){
  D_pair = data.frame(dat_pair); D_unpair = data.frame(dat_unpair)
  n_pair = nrow(D_pair); n_unpair = nrow(D_unpair); n = n_pair + n_unpair
  sign = c((dat_pair$Z1 - dat_pair$Z2 + 1)/2, dat_unpair$Z)
    
  cand_set <- rep(FALSE, n); ordered_sign <- vector(length = n); k <- 1
  while(k <= n){
    if (k %% 100 == 1) {
      pred_pair <- em_pair_sep(X1 = dat_pair$X1, Z1 = dat_pair$Z1, Y1 = dat_pair$Y1,
                               X2 = dat_pair$X2, Z2 = dat_pair$Z2, Y2 = dat_pair$Y2,
                               cand_set[1:n_pair], iter = 10)
      pred_unpair <- em_unpair_both(Y = dat_unpair$Y, Z = dat_unpair$Z, X = dat_unpair$X,
                                    cand_set = cand_set[(n_pair + 1):n], iter = 10, oracle = FALSE)
      pred = c(pred_pair, pred_unpair)
      pred = pred - min(pred) + 1 #make sure pred score strictly > 0!!!!!!
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    ordered_sign[k] <- sign[inc_ind]
    k = k + 1
  }
  k = 1:n
  S_cum = cumsum((ordered_sign - prob_assign)/sqrt(prob_assign*(1 - prob_assign)))
  m = n*cut_prop
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  }
  return(p_val)
}
  
  
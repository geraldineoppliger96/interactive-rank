
split_wilcoxon = function(dat, ub, pred_func, alternative = "two_sided", prob_assign){
  D = as.data.frame(dat)
  n = nrow(D); split_ind = sample.split(1:n, SplitRatio = 0.5)
  D = add_column(D, Z = as.factor(D$Z1 - D$Z2)); D[,c("Z1", "Z2")] = lapply(D[,c("Z1", "Z2")], factor)
  D1 = subset(D, split_ind == TRUE); D2 = subset(D, split_ind == FALSE)
  
  S = (as.numeric(D2$Z) - 1 - prob_assign)/sqrt(prob_assign*(1 - prob_assign))
  n2 = nrow(D2); k = 1
  cand_set = rep(FALSE, n2); ordered_S <- vector(length = n)
  while(k <= n2){
    if (k == 1) {
      pred = pred_func(D1, D2, cand_set)
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    ordered_S[k] <- S[inc_ind]
    k = k + 1
  }
  k = 1:n; m = n/4
  S_cum = cumsum(ordered_S)
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  }
  return(p_val)
}

combine_predict = function(D1, D2, cand_set){
  D1 = rbind(D1, subset(D2, cand_set == TRUE))
  classifier = randomForest(Z ~. - Z1 - Z2, data = D1)
  pred = predict(classifier, newdata = D2, type = "vote")[,2]
  return(pred)
}

sep_predict = function(D1, D2, cand_set){
  D1 = rbind(D1, subset(D2, cand_set == TRUE))
  
  df = data.frame(Z = as.factor(c(D1$Z1, D1$Z2)), Y = c(D1$Y1, D1$Y2),
                  X = rbind(as.matrix(select(D1, starts_with("X1."))),
                            as.matrix(select(D1, starts_with("X2.")))))
  # df1 = D1 %>% select(Z1, Y1, starts_with("X1.")); colnames(df1) = c("Z", "Y", paste("X",1:4, sep = ""))
  # df2 = D1 %>% select(Z2, Y2, starts_with("X2.")); colnames(df2) = c("Z", "Y", paste("X",1:4, sep = ""))
  # df = rbind(df1, df2)
  classifier = randomForest(Z ~ ., data = df)
  
  df = data.frame(Z = as.factor(c(D1$Z1, D1$Z2)), Y = c(D2$Y1, D2$Y2),
                  X = rbind(as.matrix(select(D2, starts_with("X1."))),
                            as.matrix(select(D2, starts_with("X2.")))))
  # df1 = D2 %>% select(Z1, Y1, starts_with("X1.")); colnames(df1) = c("Z", "Y", paste("X",1:4, sep = ""))
  # df2 = D2 %>% select(Z2, Y2, starts_with("X2.")); colnames(df2) = c("Z", "Y", paste("X",1:4, sep = ""))
  # df = rbind(df1, df2)
  pred_double = predict(classifier, newdata = df, type = "prob")[,2]
  n = nrow(D2)
  pred = pred_double[1:n] - pred_double[(n+1):(2*n)]
  pred = pred - min(pred) + 1 #predition score positive
  return(pred)
}

inverse_predict = function(D1, D2, cand_set){
  D1 = rbind(D1, subset(D2, cand_set == TRUE)); n = nrow(D2)
  df1 = D1 %>% select(Z1, Y1, starts_with("X1.")); colnames(df1) = c("Z", "Y", paste("X",1:4, sep = ""))
  df2 = D1 %>% select(Z2, Y2, starts_with("X2.")); colnames(df2) = c("Z", "Y", paste("X",1:4, sep = ""))
  df = rbind(df1, df2)
  #reg = randomForest(Y ~ ., data = df)
  qrf <- quantregForest(x = select(df, -Y), y = df$Y)
  
  df1 = D2 %>% select(Z1, Y1, starts_with("X1.")); colnames(df1) = c("Z", "Y", paste("X",1:4, sep = ""))
  df2 = D2 %>% select(Z2, Y2, starts_with("X2.")); colnames(df2) = c("Z", "Y", paste("X",1:4, sep = ""))
  df = rbind(df1, df2)
  df_sudo1 = df; df_sudo1$Z[df_sudo1$Z == 0] = 1
  df_sudo0 = df; df_sudo0$Z[df_sudo0$Z == 1] = 0
  
  #y_z0 = predict(reg, newdata = df_sudo0); y_z1 = predict(reg, newdata = df_sudo1)
  #prob_double = (df$Y - y_z0)^2 / ((df$Y - y_z1)^2  + (df$Y - y_z0)^2)
  node_sample <- predict(qrf, newdata = select(df_sudo0, -Y), what = function(t) {t}) 
  conden0 <- eden_mat(node_sample = node_sample, y_eval = df$Y)
  node_sample <- predict(qrf, newdata = select(df_sudo1, -Y), what = function(t) {t}) 
  conden1 <- eden_mat(node_sample = node_sample, y_eval = df$Y)
  pred_double = ifelse(conden0 == 0, 10, log(conden1 / conden0))
  
  pred = pred_double[1:n] - pred_double[(n + 1):(2*n)]
  pred = pred - min(pred) + 1; pred[is.na(pred)] = 1
  return(pred)
}

eden_mat = function(node_sample, y_eval){
  est_den = vector(length = length(y_eval))
  for(r in 1:nrow(node_sample)){
    hist_func = quiet(ash1(bin1(node_sample[r,],
                                ab = c(min(node_sample[r,]) - 10, max(node_sample[r,]) + 10),
                          nbin = length(node_sample[r,])/10)))
    est_den[r] = hist_func$y[which.min(abs(y_eval[r] - hist_func$x))]
  }
  return(est_den)
}


####################### unpair ########################

split_mann = function(X, Z, Y, ub, prob_assign, alternative = "two_sided"){
  D = data.frame(Z = as.factor(Z), Y = Y, X = X)
  n = nrow(D); split_ind = sample.split(1:n, SplitRatio =  0.5)
  D1 = subset(D, split_ind == TRUE); D2 = subset(D, split_ind == FALSE)
  
  n = nrow(D2); k <- 1
  cand_set <- rep(FALSE, n); ordered_Z <- vector(length = n); ordered_ind = vector(length = n)
  while(k <= n){
    if (k == 1) {
      D1 = rbind(D1, subset(D2, cand_set == TRUE))
      classifier = randomForest(Z ~., data = D1)
      pred = predict(classifier, newdata = D2, type = "vote")[,2]
      # reg = randomForest(Y ~., data = D1)
      # D2_sudo1 = D2; D2_sudo1$Z[D2_sudo1$Z == 0] = 1
      # D2_sudo0 = D2; D2_sudo0$Z[D2_sudo0$Z == 1] = 0
      # y_z0 = predict(reg, newdata = D2_sudo0, type = "prob")[,2]
      # y_z1 = predict(reg, newdata = D2_sudo1, type = "prob")[,2]
      # y = as.numeric(D2$Y)
      # pred = (y - y_z0)^2 / ((y - y_z1)^2  + (y - y_z0)^2) + 1
    }
    inc_ind <- which.max(pred*(!cand_set)); cand_set[inc_ind] <- TRUE
    ordered_Z[k] <- D2$Z[inc_ind]
    ordered_ind[k] <- inc_ind
    k = k + 1
  }
  k = 1:n; m = n/4
  #S_cum = cumsum(2*ordered_Z - 3)
  S_cum = cumsum((ordered_Z - 1 - prob_assign)/sqrt(prob_assign*(1 - prob_assign)))
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  }
  if(FALSE){
    tr_effect = (0*(D2$X.1 == 0 | D2$X.1 == 2) +
                   1*(D2$X.1 == 1))*(abs(D2$X.2)^3 - 0)*(as.numeric(D2$Z) - 1)
    tr_effect[ordered_ind[1:20]]
    sort(tr_effect, decreasing = TRUE)[1:20] # only pick up signals larger than 2
    D2$Z[ordered_ind[1:20]]
  }
  return(p_val)
}


split_mixed = function(dat_pair, dat_unpair, ub, alternative = "two_sided"){
  D_pair = data.frame(dat_pair); D_unpair = data.frame(dat_unpair)
  n_pair = nrow(D_pair)/2; n_unpair = nrow(D_unpair); n = n_pair + n_unpair

  split_ind = sample.split(1:n, SplitRatio =  0.5)
  D1 = rbind(subset(D_pair, rep(split_ind[1:n_pair], 2) == TRUE),
             subset(D_unpair, split_ind[(n_pair+1):n] == TRUE))
  D1[,c("Z")] = factor(D1[,c("Z")])
  classifier = randomForest(Z ~., data = D1)

  D2_pair = subset(D_pair, rep(split_ind[1:n_pair], 2) == FALSE)
  n_pair_test = nrow(D2_pair)/2
  pred_double = predict(classifier, newdata = D2_pair, type = "prob")[,2]
  pred_pair = pred_double[1:n_pair_test]/
    (pred_double[1:n_pair_test] + pred_double[(n_pair_test+1):(2*n_pair_test)])

  D2_unpair = subset(D_unpair, split_ind[(n_pair+1):n] == FALSE)
  pred_unpair = predict(classifier, newdata = D2_unpair, type = "prob")[,2]

  pred = c(pred_pair, pred_unpair)
  signs = c(D2_pair$Z[1:n_pair_test] - D2_pair$Z[(n_pair_test+1):(2*n_pair_test)],
            2*D2_unpair$Z - 1)
  S_cum = cumsum(signs[order(pred, decreasing = TRUE)])
  n_test = n - sum(split_ind); k = 1:n_test; m = n_test/4
  if (alternative == "two_sided") {
    p_val = 2*min(exp(-(abs(S_cum) / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "greater") {
    p_val = min(exp(-(-S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  } else if (alternative == "less") {
    p_val = min(exp(-(S_cum / (sqrt(1/2/m)*k + sqrt(m/2)))^2))
  }
  return(p_val)
}


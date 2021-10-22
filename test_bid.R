# i-bet
# dat:         a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set:    a vector of True/False to indicate whether subjects are candidates to be selected
# scale:       a factor to control the magnitude(range) of bets
# alg_type:    type of algorithm for modeling and ordering
# C_delta:     scale of treatment effect, used when alg_type = "oracle"
# binary:      whether bets are binary or continuous
# iter_round:  number of iterations we skip to update the modeling
# sum_known:   whether the total number of treatments are known
i_bid <- function(dat, cand_set,
                 scale = 0.8, alg_type = "robust", C_delta = NA, binary = FALSE,
                 iter_round = 100, sum_known = TRUE){
  n = sum(!cand_set)
  if(alg_type == "linear") {
    pred_func = em_linear
  } else if (alg_type == "robust") {
    pred_func = em_robust
  } else if (alg_type == "quadratic") {
    pred_func = em_quadratic
  } else if (alg_type == "oracle") {
    pred_func = em_oracle
  }
  
  prod_cumu <- 1; p_val <- 1; prod_seq = vector(length = 0)
  # cand_set <- rep(FALSE, n) 
  for (k in 1:n){
    if (sum_known) {
      mu_temp = sum(dat$A[!cand_set])/sum(!cand_set)
    } else {
      mu_temp = 1/2
    }
    
    if (k %% min(iter_round, n/5) == 1) {
      pred <- pred_func(dat, cand_set = cand_set, C_delta = C_delta, iter = 20)
    }
    inc_ind <- which.max(abs(pred - 1/2)*(!cand_set) + (-1)*cand_set); cand_set[inc_ind] <- TRUE
    if (binary) {
      prod_cumu = prod_cumu*
        ((2*(pred[inc_ind] > 1/2) - 1)*scale*(dat$A[inc_ind] - mu_temp) + 1)
    } else {
      prod_cumu = prod_cumu*
        ((2*(pred[inc_ind] > 1/2) - 1)*(abs(pred[inc_ind] - 1/2))*2*scale*(dat$A[inc_ind] - mu_temp) + 1)
    }
    prod_seq = c(prod_seq, prod_cumu)
    p_val <- min(p_val, 1/prod_cumu)
  }
  return(p_val)
}



# i-bet using a cross-fitting framework
# dat:         a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# scale:       a factor to control the magnitude(range) of bets
# alg_type:    type of algorithm for modeling and ordering
# C_delta:     scale of treatment effect, used when alg_type = "oracle"
# iter_round:  number of iterations we skip to update the modeling
# sum_known:   whether the total number of treatments are known
i_bid_cross <- function(dat,
                        scale = 0.8, alg_type = "robust",
                        C_delta = NA, iter_round = 100, sum_known = TRUE){
  n = length(dat$A)
  cand_set1 = rep(FALSE, n); cand_set1[sample(n, n*0.5)] = TRUE
  cand_set2 = !cand_set1
  
  p_val1 = i_bid(dat = dat, cand_set = cand_set1, scale = scale, alg_type = alg_type,
                 C_delta = C_delta, iter_round = iter_round, sum_known = sum_known)
  p_val2 = i_bid(dat = dat, cand_set = cand_set2, scale = scale, alg_type = alg_type,
                 C_delta = C_delta, iter_round = iter_round, sum_known = sum_known)
  p_val = 2*min(p_val1, p_val2)
  return(p_val)
}


# Bet and test at the end
# dat:         a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# cand_set:    a vector of True/False to indicate whether subjects are candidates to be selected
# alg_type:    type of algorithm for modeling and ordering
# scale:       a factor to control the magnitude(range) of bets
# C_delta:     scale of treatment effect, used when alg_type = "oracle"
end_bid <- function(dat, cand_set, alg_type, scale = 0.8, C_delta = NA){
  if(alg_type == "linear") {
    pred_func = em_linear
  } else if (alg_type == "robust") {
    pred_func = em_robust
  } else if (alg_type == "quadratic") {
    pred_func = em_quadratic
  } else if (alg_type == "oracle") {
    pred_func = em_oracle
  }
  
  pred <- pred_func(dat, cand_set = cand_set, C_delta = C_delta, iter = 20)
  prod_item <- ((2*(pred > 1/2) - 1)*(abs(pred - 1/2))*2*scale*(dat$A - 1/2) + 1)
  prod_cum = prod(prod_item[!cand_set])
  p_val <- 1/prod_cum
  return(p_val)
}


# Bet and test at the end, wrapped by a cross-fitting framework
# dat:         a data frame of outcome (Y), treatment assignment (A), and covariates (X.)
# alg_type:    type of algorithm for modeling and ordering
# scale:       a factor to control the magnitude(range) of bets
# C_delta:     scale of treatment effect, used when alg_type = "oracle"
end_bid_cross <- function(dat, alg_type, scale = 0.8, C_delta = NA){
  n = length(dat$A)
  cand_set1 = rep(FALSE, n); cand_set1[sample(n, n*0.5)] = TRUE
  cand_set2 = !cand_set1
  
  p_val1 = end_bid(dat = dat, cand_set = cand_set1, alg_type = alg_type, scale = scale, C_delta = C_delta)
  p_val2 = end_bid(dat = dat, cand_set = cand_set2, alg_type = alg_type, scale = scale, C_delta = C_delta)
  # p_val = 1 - pchisq(-2*(log(p_val1) + log(p_val2)), df = 4)
  p_val = 2*min(p_val1, p_val2)
  return(p_val)
}



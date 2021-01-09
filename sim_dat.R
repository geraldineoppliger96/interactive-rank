# sample_generator
# n:              total number of subjects
# m:              number of subjects with X(1) = X(2) = 1 
# treatment_type: type of treatment effect
# control_type:   type of control outcome
# C_delta:        scale of the treatment effect 
# e_D:            type of noise
# blind:          index of missing covariate
# redundant:      number of redundant variables
sample_generator = function(n, m, treatment_type, control_type, C_delta, e_D, 
                         blind = NA, redundant = NA){
  
  X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
            c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
            rnorm(n = n, sd = 1)) #a3 

  if (treatment_type == "linear") {
    Delta = C_delta*(X[,1]*X[,2] + X[,3])*2/5
  } 
  else if (treatment_type == "quadratic") {
    Delta = C_delta*(X[,3]^2 - 1)*3/5
  } 
  else if (treatment_type == "pos") {
    Delta = C_delta*(X[,1]/8 - (X[,3] > 1)/5) 
  }
  else if (treatment_type == "dense_weak") {
    Delta = C_delta*(1 - abs(sin(3*X[,3])))
  } 
  else if (treatment_type == "sparse_strong") {
    Delta = C_delta*(2*exp(X[,3])*(X[,3] > 1.5))
  }
  else if (treatment_type == "sparse_strong_smallN") {
    Delta = C_delta*(2*exp(X[,3])*(X[,3] > 1))
  }
  else if (treatment_type == "both_pos_strong") {
    Delta = C_delta*(exp(X[,3])*(X[,3] > 2) - X[,1]/2)
  }
  else if (treatment_type == "both_pos_strong_smallN") {
    Delta = C_delta*(exp(X[,3])*(X[,3] > 1) - X[,1]/2)
  }
  else if (treatment_type == "both_sparse_strong") {
    Delta = C_delta*((X[,3])^3*(abs(X[,3]) > 1))
  } 
  else if (treatment_type == "both_dense_weak") {
    Delta = C_delta*sin(3*X[,3])/5*2
  }
  
  if (control_type == "bell"){
    f = 5*rowSums(X)
  } else if (control_type == "skewed") {
    f = 2*((X[,3] < -2)*exp(-X[,3])^2)
  } 
  
  A = rbinom(n = n, size = 1, prob = 1/2)
  if (e_D == "Gaussian") {
    e = rnorm(n)
  } else if (e_D == "Beta") {
    e = rbeta(n, 2, 5)
  } else if (e_D == "Cauchy") {
    e = rcauchy(n)
  }
  Y = Delta*A + f + e
  if (!is.na(blind)) {X = X[,-blind]}
  if (!is.na(redundant)) {
    X = cbind(X, matrix(rnorm(n = redundant*n, sd = 1), ncol = redundant))
  }
  return(data.frame(Y = Y, A = A, X = X))
}



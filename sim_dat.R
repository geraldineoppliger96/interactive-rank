sample_pair_noside = function(n, p, delta, D){
  B = rbinom(n, size = 1, prob = p)
  if (D == "Gaussian") {
    Y1 = rnorm(n); Y2 = rnorm(n)
  } else if (D == "Beta") {
    Y1 = rbeta(n, 2, 5); Y2 = rbeta(n, 2, 5)
  } else if (D == "Cauchy") {
    Y1 = rcauchy(n); Y2 = rcauchy(n)
  }
  Y2 = Y2 + delta*B
  return(list(Y1 = Y1, Y2 = Y2))
}

############ pair
sample_pair_side = function(t, n, m, delta, C_delta, e_D, dat_type, blind = NA, prob){
  Z1 = rbinom(n, size = 1, prob = prob); Z2 = 1 - Z1
  
  if (e_D == "Gaussian") {
    e = rnorm(2*n)
  } else if (e_D == "Beta") {
    e = rbeta(2*n, 2, 5)
  } else if (e_D == "Cauchy") {
    e = rcauchy(2*n)
  }
  
  if (dat_type == "symmetric") {
    X1 = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
               c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
               rnorm(n = n, sd = 1)) #a3 
    X2 = cbind(X1[,c(1,2)], rnorm(n = n, sd = 1))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(X[,1]*X[,2] + X[,3])*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(X[,1]*X[,2] + X[,3])*Z + 5*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type %in% c("treatment-skewed")) {
    X1 = cbind(rbinom(n, 2, 0.9), rnorm(n))
    X2 = cbind(X1[,1], rnorm(n))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*(C_delta*Z*(0*(X[,1] == 0 | X[,1] == 2) + 1*(X[,1] == 1))*(abs(X[,2])^3 - 1)) +
          (1 - t)*(delta*Z + 1*rowSums(X))
      } else {
        Y = (C_delta*Z*(0*(X[,1] == 0 | X[,1] == 2) + 1*(X[,1] == 1))*(abs(X[,2])^3 - 1)) +
          1*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "control-skewed") {
    X1 = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
               c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
               rnorm(n = n, sd = 1), #a3
               rbeta(n = n, 2, 5)) #a4
    X2 = cbind(X1[,c(1,2)], rnorm(n = n, sd = 1), rbeta(n = n, 2, 5))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(X[,1]*X[,2] + X[,3])*Z +
          (1 - t)*(delta*Z + 2*(X[,1]*exp(X[,3])^2 + X[,4]))
      } else { 
        Y = C_delta*(X[,1]*X[,2] + X[,3])*Z + 
          2*(X[,1]*exp(X[,3])^2 + X[,4])
      }
      return(Y)
    }
  }
  else if (dat_type %in% c("new-treatment-skewed")) {
    X1 = cbind(rbinom(n, 2, 0.9), rchisq(n, df = 3))
    X2 = cbind(X1[,1], rchisq(n, df = 3))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*(C_delta*Z*(X[,2] - 1)*(X[,1] == 1)) +
          (1 - t)*(delta*Z + C_f*rowSums(X))
      } else {
        Y = (C_delta*Z*(X[,2] - 1)*(X[,1] == 1)) +
          C_f*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "nonlin_rf_extra") {
    X1 = cbind(rbinom(n, 2, 0.9), rnorm(n), rbeta(n = n, 2, 5), rnorm(n))
    X2 = cbind(X1[,1], rnorm(n),  rbeta(n = n, 2, 5), rnorm(n))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*Z*(0*(X[,1] == 0 | X[,1] == 2) + 1*(X[,1] == 1))*abs(X[,2])^3 +
          (1 - t)*(delta*Z + C_f*(X[,2] > 1))
      } else {
        Y = C_delta*Z*(0*(X[,1] == 0 | X[,1] == 2) + 1*(X[,1] == 1))*abs(X[,2])^3 +
          C_f*(X[,2] > 1)
      }
      return(Y)
    }
  }
  Y1 = out_fun(t, X1, Z1) + e[1:n]
  Y2 = out_fun(t, X2, Z2) + e[(n + 1):(2*n)]
  if (!is.na(blind)) {X1 = X1[,-blind]; X2 = X2[,-blind]}
  return(list(X1 = X1, Z1 = Z1, Y1 = Y1, X2 = X2, Z2 = Z2, Y2 = Y2))
}



##################### unpair
sample_unpair = function(t, n, m, delta, C_delta, Cf, e_D, dat_type,
                         blind = NA, redundant = NA, prob){
  Z = rbinom(n = n, size = 1, prob = prob)
  
  if (e_D == "Gaussian") {
    e = rnorm(n)
  } else if (e_D == "Beta") {
    e = rbeta(n, 2, 5)
  } else if (e_D == "Cauchy") {
    e = rcauchy(n)
  }
  
  if (dat_type == "linear_both") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(X[,1]*X[,2] + X[,3])*2/5*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  if (dat_type == "linear_both_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(X[,1]*X[,2] + X[,3])*2/5*Z + 2*((X[,3] < -2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  } 
  if (dat_type == "linear_pos") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(X[,3]^2 - 1)*3/5*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  if (dat_type == "linear_sparse") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(X[,1]*X[,2])*2/5*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  if (dat_type == "linear_varyx_big") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    X[,3] = 2*exp(X[,3])*(X[,3] > 1.5)
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*X[,3]*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  if (dat_type == "linear_varyx_small") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    X[,3] = 1 - abs(sin(3*X[,3]))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*X[,3]*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  if (dat_type == "both_big") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
               c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
               rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  else if (dat_type == "both_big_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(abs(sin(3*X[,3]))/3)*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + 2*((X[,3] < -2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  } 
  else if (dat_type == "varyx_big_wiggly") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(exp(X[,3])*(X[,3] > 2))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(exp(X[,3])*(X[,3] > 2))*Z + 5/rowSums(X)
      }
      return(Y)
    }
  } 
  else if (dat_type == "varyx_small") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(abs(sin(3*X[,3]))/3)*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(1 - abs(sin(3*X[,3])))*Z + 5*rowSums(X)
      }
      return(Y)
    }
  } 
  else if (dat_type == "varyx_small_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(abs(sin(3*X[,3]))/3)*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(1 - abs(sin(3*X[,3])))*Z + 2*((X[,3] < -2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  } 
  else if (dat_type == "varyx_big") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(exp(X[,3])*(X[,3] > 2))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(4*exp(X[,3])*(X[,3] > 2))*Z + 5*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "varyx_big_more") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(exp(X[,3])*(X[,3] > 2))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(2*exp(X[,3])*(X[,3] > 1.5))*Z + 5*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "pos_big") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(exp(X[,3])*(X[,3] > 2) - X[,1]/2)*Z + (1 - t)*(delta*Z + Cf*rowSums(X))
      } else {
        Y = C_delta*(exp(X[,3])*(X[,3] > 2) - X[,1]/2)*Z + Cf*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "both_small") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*sin(3*X[,3])*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*sin(3*X[,3])/5*2*Z + 5*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "redundant") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1), #a3 
              matrix(rnorm(n = 10*n, sd = 1), ncol = 10))
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(X[,1]*X[,2])*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(X[,1]*X[,2])*Z + 5*rowSums(X)
      }
      return(Y)
    }
  }
  else if (dat_type == "pos_big_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(X[,1]*X[,2] + X[,3])*Z +
          (1 - t)*(delta*Z + 2*(X[,1]*exp(X[,3])^2 + X[,4]))
      } else { 
        Y = C_delta*(exp(X[,3])*(X[,3] > 2) - X[,1]/2)*Z + 
          2*((X[,3] < -2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  }
  else if (dat_type == "varyx_big_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(X[,1]*X[,2] + X[,3])*Z +
          (1 - t)*(delta*Z + 2*(X[,1]*exp(X[,3])^2 + X[,4]))
      } else { 
        Y = C_delta*(4*exp(X[,3])*(X[,3] > 2))*Z + 
          2*((X[,3] < - 2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  }
  else if (dat_type == "varyx_big_more_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(exp(X[,3])*(X[,3] > 2))*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*(2*exp(X[,3])*(X[,3] > 1.5))*Z + 2*((X[,3] < - 2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  }
  else if (dat_type == "both_big_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*(X[,1]*X[,2] + X[,3])*Z +
          (1 - t)*(delta*Z + 2*(X[,1]*exp(X[,3])^2 + X[,4]))
      } else { 
        Y = C_delta*((X[,3])^3*(abs(X[,3]) > 1))*Z + 
          2*((X[,3] < -2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  }
  else if (dat_type == "both_small_control_skewed") {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
    out_fun = function(t, X, Z) {
      if (t >= 0 & t <= 1) {
        Y = t*C_delta*sin(3*X[,3])*Z + (1 - t)*(delta*Z + 5*rowSums(X))
      } else {
        Y = C_delta*sin(3*X[,3])/5*2*Z + 2*((X[,3] < -2)*exp(-X[,3])^2)
      }
      return(Y)
    }
  }
  Y = out_fun(t, X, Z) + e
  if (!is.na(blind)) {X = X[,-blind]}
  if (!is.na(redundant)) {
    X = cbind(X, matrix(rnorm(n = redundant*n, sd = 1), ncol = redundant))
  }
  return(list(X = X, Y = Y, Z = Z))
}



############## multi-sample
sample_multi_pair = function(n, p, delta, e_D){
  B = matrix(rbinom(3*n, size = 1, prob = p), ncol = 3)
  if (e_D == "Gaussian") {
    e = rnorm(3*n)
  } else if (e_D == "Beta") {
    e = rbeta(3*n, 2, 5)
  } else if (e_D == "Cauchy") {
    e = rcauchy(3*n)
  }
  
  Z_mat = t(replicate(n, sample(1:3)))
  Y_mat = B*(Z_mat == 1)*delta + B*(Z_mat == 2)*delta/2 + matrix(e, ncol = 3)

  return(list(Y = Y_mat, Z = Z_mat))
}


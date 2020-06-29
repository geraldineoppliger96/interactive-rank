ub_gua = function(alpha, m, L){
  x = sqrt(-2*m*(log(alpha))) 
  ub = x + x/2/m*(1:L - m)
  return(ub)
}

ub_ber <- function(alpha, m, L){
  u <- sqrt(-2*log(alpha)/m)
  phi <- (1 + u)/2*log(1 + u) + (1 - u)/2*log(1 - u)
  while (phi > -log(alpha)/m){
    u0 <- u
    u = u - 10^(-10)
    phi <- (1 + u)/2*log(1 + u) + (1 - u)/2*log(1 - u)
  }
  
  x = u*m
  ub <- x + (log(1/(1 - u)) - log(1 + u)) /
    (log(1/(1 - u)) + log(1 + u)) *(1:L - m)
  return(ub)
}  


dat_type = "symmetric"
n = 500
frac_pair = 0.5
prob_assign = 0.5
p = 0.05
m = 30
t = 2
delta = 0.2
C_delta = 2
Cf = 5
D = "Gaussian"
e_D = "Gaussian"
blind = NA
redundant = NA
alpha = 0.05

alg_type = "RF"
ub = ub_gua(alpha/2, n/4, n) #two-sided
ub_half = ub_gua(alpha/4, n/4, n)
R = 50
methods_pair = c("covadj-Wil", "covadj-MannW", "indep-test", "cor-test", "covadj-Wil-value",
                 "split-Wil-sep", "split-Wil-combine", 
                 "pseudo-Wil-combine-rank", "pseudo-Wil-sep-rank", 
                 "pseudo-Wil-combine-value", "pseudo-Wil-sep-value")
methods_unpair = c("covadj-Wil", "covadj-Wil-value",
                   "covadj-MannW", "indep-test", "cor-test",
                   "split-MannW", "interactive-default",
                   "rank-pseudo_res", "value-pseudo_res",
                   "rank-pseudo_outcome",
                   "rank-signed-pseudo1_res", "value-signed-pseudo1_res",
                   "rank-signed-pseudo2_res", "value-signed-pseudo2_res",
                   "rank-signed-pseudo3_res", "value-signed-pseudo3_res",
                   "rank-signed-pseudo3_outcome",
                   "rank-mask_res", "value-mask_res",
                   "rank-false_res", "value-false_res", "mann-false_res",
                   "rank-false_outcome",
                   "rank-aver_res", "value-aver_res", "mann-aver_res",
                   "rank-diff_res", "value-diff_res")
n_permute = 200

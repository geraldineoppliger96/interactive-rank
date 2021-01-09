alpha = 0.05
n = 500
m = 30
e_D = "Gaussian"
blind = NA
redundant = NA
R = 500
iter_round = 100

methods_interactive = c("CovAdj-Wilcoxon-linear", "CovAdj-Wilcoxon-robust", "CovAdj-Wilcoxon-quadratic",
                        "linear-CATE-test",
                        "i-Wilcoxon-linear", "i-Wilcoxon-robust",
                        "i-Wilcoxon-quadratic")

n_permute = 200
alg_type = "RF"
methods_var_Wilcoxon =  c("R(X)",
                          "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
                          "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
                          "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")

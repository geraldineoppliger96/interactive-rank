alpha = 0.05
n = 500
m = 30
e_D = "Gaussian"
blind = NA
redundant = NA
R = 500
iter_round = 100
scale = 0.8
split_pct = 0.1

methods_interactive = c("CovAdj-Wilcoxon-linear", "CovAdj-Wilcoxon-robust",
                        "CovAdj-Wilcoxon-quadratic", "linear-CATE-test",
                        "i-Wilcoxon-robust-signedA_adapt",
                        "i-Wilcoxon-quadratic-signedA_adapt",
                        "i-bid-linear", "i-bid-robust", "i-bid-quadratic",
                        "i-bid-linear-bi", "i-bid-robust-bi", "i-bid-quadratic-bi",
                        "i-bid-cross-linear", "i-bid-cross-robust", "i-bid-cross-quadratic")

n_permute = 200
alg_type = "RF"
methods_var_Wilcoxon =  c("R(X)",
                          "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
                          "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
                          "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")

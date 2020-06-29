source("setup.R")
source("input.R", local = TRUE)
args <- commandArgs(trailingOnly = TRUE)

mode = args[1]
for(single_mode in strsplit(mode, split = ";")[[1]]){
  para_vary = strsplit(single_mode, split = ":")[[1]]
  if (para_vary[1] %in% c("dat_type", "e_D")){
    assign(para_vary[1], para_vary[2])
  } else {
    assign(para_vary[1], as.numeric(para_vary[2]))
  }
}
C_delta = as.integer(args[2])
index = as.integer(args[3])
path = args[4]

rejection = vector(length = length(methods_unpair)); names(rejection) = methods_unpair

dat = sample_unpair(t = t, n = n, m = m, delta = delta, prob = prob_assign, redundant = redundant,
                    C_delta = C_delta, Cf = Cf, e_D = e_D, dat_type = dat_type)
df = data.frame(Y = dat$Y, dat$X);
if (alg_type %in% c("linear")) {
  e = lm(Y ~ .^2 - 1, data = df)$residuals
} else if(alg_type == "RF"){
  e = dat$Y - randomForest(Y ~ ., data = df)$predicted
} 
if ("MannW" %in% methods_unpair) {
  rejection["MannW"] = wilcox.test(dat$Y ~ dat$Z, paired = FALSE)$p.value 
}
if ("Wil" %in% methods_unpair) {
  rejection["Wil"] =
    wilcox.test(dat$Y*(2*dat$Z - 1), rep(0, n), paired = TRUE)$p.value 
}
if ("Wil-value" %in% methods_unpair) {
  rejection["Wil-value"] = 
    permute_wilcoxon(Y_diff = dat$Y, Z_diff = (2*dat$Z - 1), n_permute = n_permute) 
}
if ("covadj-MannW" %in% methods_unpair) {
  rejection["covadj-MannW"] = wilcox.test(e ~ dat$Z, paired = FALSE, alternative = "less")$p.value 
}
if ("covadj-Wil" %in% methods_unpair) {
  rejection["covadj-Wil"] = 
    wilcox.test((2*dat$Z - 1)*e, rep(0, n), paired = TRUE, alternative = "greater")$p.value 
}
if ("covadj-Wil-value" %in% methods_unpair) {
  rejection["covadj-Wil-value"] = 
    permute_wilcoxon(Y_diff = e, Z_diff = (2*dat$Z - 1), n_permute = n_permute) 
}
if ("CATE-test" %in% methods_unpair) {
  rejection["CATE-test"] = CATE_test(X = dat$X, Z = dat$Z, Y = dat$Y, alpha = alpha)
}
if ("indep-test" %in% methods_unpair) {
  rejection["indep-test"] = dcov.test(e, dat$Z, R = 100)$p.value 
}
if ("cor-test" %in% methods_unpair) {
  rejection["cor-test"] = cor.test(e, dat$Z)$p.value
}
if ("split-MannW" %in% methods_unpair) {
  rejection["split-MannW"] =
    split_mann(X = dat$X, Z = dat$Z, Y = dat$Y, ub = ub, prob_assign = prob_assign) 
}
if ("interactive-default" %in% methods_unpair) {
  rejection["interactive-default"] = 
    i_mann_pval(X = dat$X, Z = dat$Z, Y = dat$Y, prob_assign = prob_assign, type = "default")
}
if (any(grepl("res", methods_unpair))) {
  df = data.frame(Y = e, Z = dat$Z, X = dat$X)
  p_res_vector = batch_pseudo_wilcoxon(dat = df, alg_type = alg_type, n_permute = n_permute)
  for(avail_method in names(p_res_vector)){
    if (paste(avail_method, "_res", sep = "") %in% methods_unpair) {
      rejection[paste(avail_method, "_res", sep = "")] =
        p_res_vector[avail_method]
    }
  }
}
if (any(grepl("outcome", methods_unpair))) {
  df = data.frame(Y = dat$Y, Z = dat$Z, X = dat$X)
  p_outcome_vector = batch_pseudo_wilcoxon(df, alg_type = alg_type, n_permute = n_permute)
  for(avail_method in names(p_outcome_vector)){
    if (paste(avail_method, "_outcome", sep = "") %in% methods_unpair) {
      rejection[paste(avail_method, "_outcome", sep = "")] =
        p_outcome_vector[avail_method]
    }
  }
}

write.table(data.frame(rejection), file = sprintf("%s/%d_%d.txt", path, C_delta, index), 
            row.names = FALSE, col.names = FALSE)


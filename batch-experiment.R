source("setup.R")
dir.create("result")

####################################################################
######################### i-Wilcoxon test ##########################
####################################################################

################### Figure 2
dat_type = "linear_both"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "methods_unpair",
                         value = c("covadj-Wil", "covadj-MannW", "CATE-test", 
                                   "interactive",
                                   "interactive-default",
                                   "interactive-naive")),
                   list(name = "dat_type", value = dat_type),
                   list(name = "alg_type", value = "linear"),
                   list(name = "n_permute", value = 100),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "methods_unpair",
                        value = c("covadj-Wil", "covadj-MannW", "CATE-test", 
                                  "interactive-default",
                                  "interactive-naive")),
                   list(name = "dat_type", value = dat_type),
                   list(name = "alg_type", value = "linear"),
                   list(name = "blind", value = 3),
                   list(name = "n_permute", value = 100),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/blind_", dat_type, ".Rdata", sep = ""))

################### Figure 3
dat_type = "linear_both_control_skewed"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "methods_unpair",
                        value = c("covadj-Wil", "covadj-MannW", "CATE-test", 
                                  "interactive",
                                  "interactive-default",
                                  "interactive-naive")),
                   list(name = "dat_type", value = dat_type),
                   list(name = "alg_type", value = "linear"),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))


################### Figure 4
dat_type = "linear_pos"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "methods_unpair",
                        value = c("covadj-Wil", "covadj-MannW", "CATE-test", 
                                  "interactive",
                                  "interactive-default",
                                  "interactive-naive")),
                   list(name = "dat_type", value = dat_type),
                   list(name = "alg_type", value = "linear"),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))



####################################################################
######################### variants of Wilcoxon #####################
####################################################################

############## Figure 5
dat_type = "varyx_small"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))

dat_type = "varyx_small"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/cauchy_", dat_type, ".Rdata", sep = ""))

dat_type = "varyx_big_more"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))

dat_type = "varyx_small_control_skewed"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))


############## Figure 6
dat_type = "varyx_big_more_control_skewed"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))


############## Figure 7
dat_type = "varyx_big_more"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/cauchy_", dat_type, ".Rdata", sep = ""))



############## Figure 7
dat_type = "pos_big"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))

dat_type = "both_big"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))

dat_type = "both_small"; max_sig = 5
Cd_seq = seq(0, max_sig, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  print(Cd)
  para_vary = list(list(name = "dat_type", value = dat_type),
                   list(name = "C_delta", value = Cd),
                   list(name = "R", value = 500))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste("result/", dat_type, ".Rdata", sep = ""))




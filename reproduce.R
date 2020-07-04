source("setup.R")
if (!file.exists("result")) {
  dir.create("result")
}
############################ Figure 2 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,".Rdata", sep = ""))

treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "blind", value = 3),
                   list(name = "methods_interactive",
                        value = c("CovAdj-Wilcoxon-linear", "CovAdj-Wilcoxon-robust", 
                                  "linear-CATE-test",
                                  "i-Wilcoxon-linear", "i-Wilcoxon-robust")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_blind.Rdata", sep = ""))


############################ Figure 3 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_control_skewed.Rdata", sep = ""))

############################ Figure 4 ##################################
treatment_type = "quadratic"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,".Rdata", sep = ""))

############################ Figure 10 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"skewed_cauchy.Rdata", sep = ""))


########################################################################
######################### variant Wilcoxon #############################
########################################################################

############# dense-weak effect
###############################
treatment_type = "dense_weak"; max_scale = 5
Cd_seq = seq(0, max_scale, length.out = 6)

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, ".Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_control_skewed.Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_cauchy.Rdata", sep = ""))


########## sparse-strong effect
###############################
treatment_type = "sparse_strong"; max_scale = 5
Cd_seq = seq(0, max_scale, length.out = 6)

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, ".Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_control_skewed.Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_cauchy.Rdata", sep = ""))


############## two-sided effect
###############################
treatment_type = "both_pos_strong"; max_scale = 5
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, ".Rdata", sep = ""))

treatment_type = "both_sparse_strong"; max_scale = 5
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, ".Rdata", sep = ""))

treatment_type = "both_dense_weak"; max_scale = 5
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, ".Rdata", sep = ""))




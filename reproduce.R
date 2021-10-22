source("setup.R")
if (!file.exists("result")) {
  dir.create("result")
}


methods_updated = c("CovAdj-Wilcoxon-linear", "CovAdj-Wilcoxon-robust",
                    "CovAdj-Wilcoxon-quadratic", "linear-CATE-test",
                    "i-Wilcoxon-robust-signedA_adapt",
                    "i-Wilcoxon-quadratic-signedA_adapt",
                    "i-bid-linear", "i-bid-robust", "i-bid-quadratic",
                    "i-bid-linear-bi", "i-bid-robust-bi", "i-bid-quadratic-bi",
                    "i-bid-cross-linear", "i-bid-cross-robust", "i-bid-cross-quadratic")

############################ Figure 2 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/i_bet_", treatment_type,"_bell.Rdata", sep = ""))



############################ Figure 3 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/i_bet_", treatment_type,"_skewed.Rdata", sep = ""))




############################ Figure 4 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/i_bet_", treatment_type,"_cauchy.Rdata", sep = ""))


############################ Figure 5 ##################################
treatment_type = "quadratic"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/i_bet_", treatment_type,"_bell.Rdata", sep = ""))


if(0){
############################ Figure 12 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "methods_interactive",
                        value = c("i-Wilcoxon-robust-signedA", 
                                  "i-Wilcoxon-robust-signedA_adapt")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_adapt.Rdata", sep = ""))


treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "methods_interactive",
                        value = c("i-Wilcoxon-robust-signedA", 
                                  "i-Wilcoxon-robust-signedA_adapt")),
                   list(name = "iter_round", value = 20),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"control_skewed_adapt.Rdata", sep = ""))

treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "methods_interactive",
                        value = c("i-Wilcoxon-robust-signedA", 
                                  "i-Wilcoxon-robust-signedA_adapt")),
                   list(name = "iter_round", value = 20),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"cauchy_adapt.Rdata", sep = ""))

treatment_type = "quadratic"
max_scale = 5/3; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "methods_interactive",
                        value = c("i-Wilcoxon-quadratic-signedA", 
                                  "i-Wilcoxon-quadratic-signedA_adapt")),
                   list(name = "iter_round", value = 20),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_adapt.Rdata", sep = ""))

############################ Figure 14 ##################################
treatment_type = "linear"
max_scale = 5; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_cauchy_updated.Rdata", sep = ""))


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



#########################################################################
########################### small sample size  ##########################
#########################################################################


#########################################################################
########################### Figure 15 ###################################
#########################################################################


############################ Figure 15 (a) ##################################
treatment_type = "linear"
max_scale = 10; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_control_skewed_updated_small.Rdata", sep = ""))


############################ Figure 15 (b) ##################################
treatment_type = "linear"
max_scale = 10; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_cauchy_updated_small.Rdata", sep = ""))


############################ Figure 15 (c) ##################################
treatment_type = "quadratic"
max_scale = 10; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "methods_interactive",
                        value = methods_updated),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_updated_small.Rdata", sep = ""))



############################ Figure 16 ##################################
treatment_type = "linear"
max_scale = 10; Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "methods_interactive",
                        value = c("i-Wilcoxon-robust-signedA", "i-Wilcoxon-oracle")),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_interactive(para_vary)
}
save(result, file=paste("result/", treatment_type,"_cauchy_oracle_small.Rdata", sep = ""))





########################################################################
######################### Figure 17 ####################################
########################################################################

############# dense-weak effect
###############################
treatment_type = "dense_weak"; max_scale = 10
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_small.Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_control_skewed_small.Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_cauchy_small.Rdata", sep = ""))


########## sparse-strong effect
###############################
treatment_type = "sparse_strong_smallN"; max_scale = 10
Cd_seq = seq(0, max_scale, length.out = 6)

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_small.Rdata", sep = ""))

result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "skewed"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_control_skewed_small.Rdata", sep = ""))


treatment_type = "sparse_strong_smallN"; max_scale = 10
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "e_D", value = "Cauchy"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_cauchy_small.Rdata", sep = ""))


############## two-sided effect
###############################
treatment_type = "both_pos_strong_smallN"; max_scale = 10
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_small.Rdata", sep = ""))


treatment_type = "both_sparse_strong"; max_scale = 10
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_small.Rdata", sep = ""))


treatment_type = "both_dense_weak"; max_scale = 10
Cd_seq = seq(0, max_scale, length.out = 6)
result = list()
for (Cd in Cd_seq) {
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "control_type", value = "bell"),
                   list(name = "n", value = 50), 
                   list(name = "m", value = 10),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_var_Wilcoxon(para_vary)
}
save(result, file=paste("result/", treatment_type, "_small.Rdata", sep = ""))
}
}




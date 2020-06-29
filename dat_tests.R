source("setup.R")
load("dat3/pphdata.RData")
all_dat = subset(table, table$study == 2); n = nrow(all_dat)
covariates = all_dat[,c("age", "q7_1", "q6_1", "par_cat", "q13_1",
             "prehb_a", "q5_2", "q5a_2", "anysuturing", "q8_2", "q9_2", "q10_2",
             "timeplacenta", "q2_3")]
names(covariates) = c("age", "marital", "edu", "parity",
               "dur_gestation", "pre_Hb", "lab_ind", "lab_aug", 
               "anysuturing", "cord_clamp","cord_traction","ult_massage",
               "timeplacenta", "bldloss_PPH")
col_factors = c("marital", "edu", "parity", "lab_ind", "lab_aug", 
                "anysuturing", "cord_clamp","cord_traction","ult_massage")
covariates[col_factors] = lapply(covariates[col_factors], factor)

exp_per_outcome = function(name_out, covariates, all_dat, res_type, methods){
  #if(name_out %in% c("posthb_a", "activebldcontrol")) {alter = "less"} else {alter = "greater"}
  dat = cbind(all_dat[,c("txtarm",name_out)], covariates); names(dat)[1:2] = c("Z", "Y")
  dat =  drop_na(dat); n = nrow(dat); dat$Z = dat$Z - 1
  
  if (res_type %in% c("linear")) {
    e = lm(Y ~ .^2 - 1 - Z, data = dat)$residuals
  } else if(res_type == "nonlinear"){
    e = dat$Y - randomForest(Y ~ . - Z, data = dat)$predicted
  } 
  
  p_vec = vector(length = length(methods)); names(p_vec) = methods
  if ("MannW" %in% methods) {
    p_vec["MannW"] = wilcox.test(dat$Y ~ dat$Z, paired = FALSE)$p.value
  }
  if ("Wil" %in% methods) {
    p_vec["Wil"] =
      wilcox.test(dat$Y*(2*dat$Z - 1), rep(0, n), paired = TRUE)$p.value
  }
  if ("covadj-MannW" %in% methods) {
    p_vec["covadj-MannW"] = wilcox.test(e ~ dat$Z, paired = FALSE)$p.value
  }
  if ("covadj-Wil" %in% methods) {
    p_vec["covadj-Wil"] = 
      wilcox.test((2*dat$Z - 1)*e, rep(0, n), paired = TRUE)$p.value
  }
  if ("CATE-test" %in% methods) {
    p_vec["CATE-test"] = CATE_test(X = dat$X, Z = dat$Z, Y = dat$Y, alpha = alpha)
  }
  if ("indep-test" %in% methods) {
    p_vec["indep-test"] = dcov.test(e, dat$Z, R = 100)$p.value
  }
  if ("cor-test" %in% methods) {
    p_vec["cor-test"] = cor.test(e, dat$Z)$p.value
  }
  if ("split-MannW" %in% methods) {
    p_vec["split-MannW"] =
      split_mann(X = dat[,-c(1,2)], Z = dat$Z, Y = dat$Y, ub = ub)
  }
  if ("sudo-Wil-rank" %in% methods) {
    p_vec["sudo-Wil-rank"] =
      sudo_wilcoxon(dat, n_permute = 100, 
                    rank_type = TRUE, alg_type = "sudo", max_nodes = NULL)
  }
  if ("sudo-Wil-value" %in% methods) {
    p_vec["sudo-Wil-value"] =
      sudo_wilcoxon(dat, n_permute = 100, 
                    rank_type = FALSE, alg_type = "sudo", max_nodes = NULL)
  }
  return(p_vec)
}



#continuous
outcomes = c("timebldstop", "addbldloss", "overallbld", "posthb_a")
methods = c("Wil", "MannW", "covadj-Wil", "covadj-MannW", "indep-test", "cor-test",
            "split-MannW", "sudo-Wil-rank", "sudo-Wil-value")
prop_seq = seq(0.6, 1, 0.1)
R = 50
result = list()
for (prop in prop_seq) {
  print(prop)
  wrapper_func = function(i) {
    print(i)
    p_mat = matrix(nrow = 0, ncol = length(methods))
    select_ind = sample(n, n*0.5)
    #rownames(p_mat) = outcomes; 
    colnames(p_mat) = methods
    for(name_out in outcomes){
      temp = exp_per_outcome(name_out, covariates[select_ind,], all_dat[select_ind,],
                             res_type = "nonlinear", methods = methods)
      p_mat = rbind(temp,p_mat); rownames(p_mat)[1] = name_out
    }
    return(p_mat)
  }
  result[[as.character(prop)]] = lapply(1:R, wrapper_func)
}
save(result, file="result/realdat_varysize.Rdata")


mean_power = sapply(result, function(x){colMeans(matrix(unlist(x), byrow = TRUE, nrow = length(x)))})
mean_sd = sapply(result, function(x){colSes(matrix(unlist(x), byrow = TRUE, nrow = length(x)))})
df_power = data.frame(methods = rep(rep(colnames(result[[1]][[1]]), each = 4), 3),
                      outcomes = rep(rownames(result[[1]][[1]]), 27),
                      prop = rep(prop_seq, each = 36),
                      pval = as.vector(mean_power))
mode = "timebldstop"
p = ggplot(data = subset(df_power, outcomes == mode &
                           !(methods %in% c("indep-test","covadj-Wil","MannW"))),
           aes(x = prop, y = pval, group = methods, fill = methods)) +
  geom_line(aes(linetype = methods, color = methods), size = 0.8) +
  geom_point(aes(shape = methods, color = methods), size = 2.5) +
  # scale_linetype_manual(values = c("dotdash", "dashed", "twodash", "solid")) +
  # scale_shape_manual(values = c(24, 22, 25, 21)) +
  # scale_color_manual(values = c( "#00FF00","#FF3366", "#FF9900",  "#0099FF")) +
  # scale_fill_manual(values = c( "#00FF00",  "#FF3366", "#FF9900",  "#0099FF")) +
  theme(legend.title = element_blank(), #legend.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "bottom", legend.text = element_text(size = 8)) +
  xlab("propotion of sample") + ylab("p-values") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
plot(p)
ggsave(filename = paste("figure/realpval_",mode,".png", sep = ""), plot = p, width = 4, height = 3.6)

leg <- get_legend(p); p2 = as_ggplot(leg)
plot(p2)
ggsave(filename = "figure/legend_real.png", plot = p2, width = 5, height = 0.5)

# #binary
# outcomes = c("activebldcontrol", "add300loss", "add500loss", "add1000loss", "hbdrop2", "hbdrop3",
#              "anyuterotonics", "anytransfus", "anyhysterectomy", "anyexploration", "anybimanual", "anyothrsurgery")
# p_mat_bi = matrix(nrow = length(outcomes), ncol = 6)
# rownames(p_mat_bi) = outcomes; colnames(p_mat_bi) = c("classical", "covadj", "cate", "interactive", "split", "number_selected")
# for(name_out in outcomes){
#   p_mat_bi[name_out,] = exp_per_outcome(name_out, covariates, all_dat, type = "bi")
#   print(p_mat_bi[name_out,])
# }
# 
# outcome = "overallbld" #"timebldstop", "addbldloss", "overallbld", "posthb_a",
# dat = dat[,c("txtarm", "age", "par_cat", "prehb_a", outcome)]
# names(dat)[5] = "outcome"
# dat =  drop_na(dat)
# 
# 
# 
# #mann-whitney two-sided
# wilcox.test(outcome ~ txtarm, paired = FALSE, data = dat)$p.value #bldloss
# t.test(outcome ~ txtarm, data = dat)$p.value                      #timebldstop, posthb_a
# #one-sided
# wilcox.test(outcome ~ txtarm, alternative = "greater", paired = FALSE, data = dat)$p.value 
# t.test(outcome ~ txtarm, data = dat, alternative = "less")$p.value
# order_mann_pval(Z = dat$txtarm - 1, Y = dat$outcome, alternative = "greater")
# 
# #cov-adjusted mann-whitney
# e = lm(outcome ~ .^2 - 1, data = dat[,-c(1,2,3)])$residuals
# wilcox.test(e ~ dat$txtarm, paired = FALSE, alternative = "less")$p.value 
# 
# #CATE-test
# CATE_test_pval(Z = dat$txtarm - 1, Y = dat$outcome,
#                X = data.frame(dat$age, dat$par_cat))
# 
# #i-mann-whitney
# summary(lm(outcome ~ .^2, data = dat[,-c(1)]))
# select_model = randomForest(outcome ~., data = dat[,-c(1)]); varImpPlot(select_model)
# select_cov = dat[,c("bldloss_PPH", "pre_Hb", "age", "timeplacenta")]
# select_cov = dat[,c("pre_Hb", "age")]
# i_mann_pval(Z = dat$treatment - 1, Y = dat$outcome,
#             X = select_cov, alternative = "greater") 
# i_mann_pval(Z = dat$txtarm - 1, Y = dat$outcome,
#             X = data.frame(dat$age, dat$par_cat), alternative = "greater") 
# i_mann_pval(Z = dat$txtarm - 1, Y = dat$outcome, X = data.frame(dat$prehb_a), alternative = "less")
# #dat$iv_fluids, dat$prehb_a
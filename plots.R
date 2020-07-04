if (!file.exists("figure")) {
  dir.create("figure")
}
source("input.R")
library(ggpubr)
colSes = function(m) {apply(m, 2, sd)/sqrt(nrow(m))} #standard error
plot_df = function(mean_power, compare_methods, legend_name){
  df_power = data.frame(mu_seq = rep(as.numeric(colnames(mean_power)),
                                     each = nrow(mean_power)),
                        power = as.vector(mean_power),
                        grp = rep(legend_name,
                                  ncol(mean_power)))
  p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
             aes(x = mu_seq, y = power, group = grp, fill = grp)) +
    geom_line(aes(linetype = grp, color = grp), size = 0.8) +
    geom_hline(yintercept=0.05) +
    geom_point(aes(shape = grp, color = grp), size = 2.5) +
    theme(legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = c(0.8,0.6), legend.text = element_text(size = 6)) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    xlab("Scale of treatment effect") + ylab("power") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
  plot(p) 
  return(p)
}

############## Figure 2
mode = "linear"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[c(1,5)] = c("CovAdj-Wilcoxon", "i-Wilcoxon")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


mode = "linear_blind"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[c(1,5)] = c("CovAdj-Wilcoxon", "i-Wilcoxon")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

############## Figure 3
mode = "linear_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[5] = c("i-Wilcoxon-original")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon-robust", "linear-CATE-test",
                    "i-Wilcoxon-original", "i-Wilcoxon-robust")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

############## Figure 4
mode = "quadratic"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]])
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon-quadratic", "linear-CATE-test",
                    "i-Wilcoxon-robust", "i-Wilcoxon-quadratic")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

############## Figure 10
mode = "linear_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[5] = c("i-Wilcoxon-original")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon-robust", "linear-CATE-test",
                    "i-Wilcoxon-original", "i-Wilcoxon-robust")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)



############### Figure 5(a)
mode = "dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
mode = "dense_weak_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power_cauchy = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R(x)-Gaussian", "R(X)-Cauchy", "R(X, 1-A)-Gaussian", "R(X, 1-A)-Cauchy")

legend_name = factor(compare_methods, levels = compare_methods)
df_combine = data.frame(mu_seq = rep(0:5, each = 4),
                     power = as.vector(rbind(mean_power[c(1,2),],
                                             mean_power_cauchy[c(1,2),])[c(1,3,2,4),]),
                     grp = rep(legend_name,
                               ncol(mean_power)))
p = ggplot(data = df_combine,
           aes(x = mu_seq, y = power, group = grp, fill = grp)) +
  geom_line(aes(linetype = grp, color = grp), size = 0.8) +
  geom_hline(yintercept=0.05) +
  geom_point(aes(shape = grp, color = grp), size = 2.5) +
  scale_shape_manual(values = c(21, 21, 24, 24)) +
  scale_color_manual(values = c("red1", "pink1", "green3", "yellow4")) +
  scale_fill_manual(values = c("red1", "pink1", "green3", "yellow4")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = c(0.77, 0.21), legend.text = element_text(size = 8)) +
  xlab("Scale of treatment effect") + ylab("power") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
plot(p)
ggsave(filename = paste("figure/outcome_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


legend_name = factor(methods_var_Wilcoxon, levels = methods_var_Wilcoxon)
############### Figure 5(b)
mode = "sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R(X)", "R(X, 1-A)")

df_power = data.frame(mu_seq = rep(0:5, each = nrow(mean_power)),
                      power = as.vector(mean_power),
                      grp = rep(legend_name,
                                ncol(mean_power)))
p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
           aes(x = mu_seq, y = power, group = grp, fill = grp)) +
  geom_line(aes(linetype = grp, color = grp), size = 0.8) +
  geom_hline(yintercept=0.05) +
  geom_point(aes(shape = grp, color = grp), size = 2.5) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = c("red1", "green3")) +
  scale_fill_manual(values = c("red1", "green3")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = c(0.16, 0.86), legend.text = element_text(size = 8)) +
  xlab("Scale of treatment effect") + ylab("power") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
plot(p)
ggsave(filename = paste("figure/outcome_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


############### Figure 5(c)
mode = "dense_weak_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R(X)", "R(X, 1-A)")

df_power = data.frame(mu_seq = rep(0:5, each = nrow(mean_power)),
                      power = as.vector(mean_power),
                      grp = rep(legend_name,
                                ncol(mean_power)))
p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
           aes(x = mu_seq, y = power, group = grp, fill = grp)) +
  geom_line(aes(linetype = grp, color = grp), size = 0.8) +
  geom_hline(yintercept=0.05) +
  geom_point(aes(shape = grp, color = grp), size = 2.5) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = c("red1", "green3")) +
  scale_fill_manual(values = c("red1", "green3")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = c(0.16, 0.86), legend.text = element_text(size = 8)) +
  xlab("Scale of treatment effect") + ylab("power") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
plot(p)
ggsave(filename = paste("figure/outcome_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


############### Figure 6(a)
mode = "dense_weak_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R(X, 1-A)",  "R - hat(R)(X, 1 - A)")
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))

df_power = data.frame(mu_seq = rep(0:5, each = nrow(mean_power)),
                      power = as.vector(mean_power),
                      grp = rep(legend_name,
                                ncol(mean_power)))
p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
           aes(x = mu_seq, y = power, group = grp, fill = grp)) +
  geom_line(aes(linetype = grp, color = grp), size = 0.8) +
  geom_hline(yintercept=0.05) +
  geom_point(aes(shape = grp, color = grp), size = 2.5) +
  scale_shape_manual(values = c(21, 24), labels = names) +
  scale_color_manual(values = c("green3", "green4"), labels = names) +
  scale_fill_manual(values = c("green3", "green4"), labels = names) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = names) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = c(0.2, 0.85), legend.text = element_text(size = 8)) +
  xlab("Scale of treatment effect") + ylab("power") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
plot(p)
ggsave(filename = paste("figure/outcome_compare_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


############### Figure 6(b)
mode = "sparse_strong_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R(X, 1-A)",  "R - hat(R)(X, 1 - A)")
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))

df_power = data.frame(mu_seq = rep(0:5, each = nrow(mean_power)),
                      power = as.vector(mean_power),
                      grp = rep(legend_name,
                                ncol(mean_power)))
p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
           aes(x = mu_seq, y = power, group = grp, fill = grp)) +
  geom_line(aes(linetype = grp, color = grp), size = 0.8) +
  geom_hline(yintercept=0.05) +
  geom_point(aes(shape = grp, color = grp), size = 2.5) +
  scale_shape_manual(values = c(21, 24), labels = names) +
  scale_color_manual(values = c("green3", "green4"), labels = names) +
  scale_fill_manual(values = c("green3", "green4"), labels = names) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = names) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = c(0.2, 0.85), legend.text = element_text(size = 8)) +
  xlab("Scale of treatment effect") + ylab("power") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
plot(p)
ggsave(filename = paste("figure/outcome_compare_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)



plot_df = function(mean_power, compare_methods, label_names) {
  legend_name = factor(methods_var_Wilcoxon, levels = methods_var_Wilcoxon)
  df_power = data.frame(mu_seq = rep(0:5, each = nrow(mean_power)),
                        power = as.vector(mean_power),
                        grp = rep(legend_name,
                                  ncol(mean_power)))
  p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
             aes(x = mu_seq, y = power, group = grp, fill = grp)) +
    geom_line(aes(linetype = grp, color = grp), size = 0.8) +
    geom_hline(yintercept=0.05) +
    geom_point(aes(shape = grp, color = grp), size = 2.5) +
    scale_shape_manual(values = c(21, 24, 22, 3), labels = label_names) +
    scale_color_manual(values = c("red1", "green4", "deepskyblue2", "purple"), labels = label_names) +
    scale_fill_manual(values = c("red1", "green4", "deepskyblue2", "purple"), labels = label_names) +
    scale_linetype_manual(values = c("solid", "22", "42", "44"), labels = label_names) +
    theme(legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = "none", legend.text = element_text(size = 8)) +
    xlab("Scale of treatment effect") + ylab("power") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) 
  plot(p)
  return(p)
}
################ Figure 7
compare_methods = c("R(X)",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R(X), R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))

mode = "sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


################### Figure 8
compare_methods = c("R(X)",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R(X), R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))

mode = "dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


################### Figure 9
compare_methods = c("R(X)",  "R - hat(R)(X, 1 - A)",
                    "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
                    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
label_names = expression(R(X), R - hat(R)(X, 1 - A), 
                         paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|"),
                         S %.% (paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|")))

mode = "both_pos_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)







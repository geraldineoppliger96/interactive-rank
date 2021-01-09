if (!file.exists("figure")) {
  dir.create("figure")
}
source("setup.R")
source("input.R")
library(ggpubr)
colSes = function(m) {apply(m, 2, sd)/sqrt(nrow(m))} #standard error
plot_df = function(mean_power, compare_methods,
                   legend_name, pos = c(0.75, 0.6), legend_size = 8){
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
          legend.position = pos, legend.text = element_text(size = legend_size)) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    xlab("Scale of treatment effect") + ylab("power") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
  plot(p) 
  return(p)
}


############## Figure 2
mode = "linear_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[c(1,8)] = c("CovAdj-Wilcoxon", "i-Wilcoxon")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

############## Figure 3
mode = "linear_control_skewed_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[c(6,8,9)] = c("i-Wilcoxon-robust_old", "i-Wilcoxon-original", "i-Wilcoxon-robust")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon-robust", "linear-CATE-test",
                    "i-Wilcoxon-original", "i-Wilcoxon-robust")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = c(0.7, 0.6))
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

############## Figure 4
mode = "quadratic_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); 
legend_name[c(6, 7, 9, 10)] = c("i-Wilcoxon-robust_old", "i-Wilcoxon-quadratic_old",
                                "i-Wilcoxon-robust", "i-Wilcoxon-quadratic")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon-quadratic", "linear-CATE-test",
                    "i-Wilcoxon-robust", "i-Wilcoxon-quadratic")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = c(0.75, 0.5), legend_size = 6)
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
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
ggsave(filename = paste("figure/outcome_",mode,".eps", sep = ""),
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
ggsave(filename = paste("figure/outcome_",mode,".eps", sep = ""),
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
ggsave(filename = paste("figure/outcome_",mode,".eps", sep = ""),
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
ggsave(filename = paste("figure/outcome_compare_",mode,".eps", sep = ""),
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
ggsave(filename = paste("figure/outcome_compare_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)



plot_df_var = function(mean_power, compare_methods, label_names) {
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
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)


################### Figure 8
compare_methods = c("R(X)",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R(X), R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))

mode = "dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".eps", sep = ""),
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
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_var(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)




#########################################################################
########################### Appendix ###################################
#########################################################################

########################### Figure 12 ###################################
mode = "linear_adapt"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "linear_control_skewed_adapt"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "linear_cauchy_adapt"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "quadratic_adapt"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original weights)", "i-Wilcoxon (new weights)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)


#########################################################################
########################### Figure 13 ###################################
#########################################################################
legend_name = c("R(X)", "Bonferroni", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|", "Geometric",
                "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)", "Harmonic")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = legend_name
label_names = expression(R(X), "Bonferroni",
                         paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|"),
                         "Geometric",
                         S %.% (paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|")),
                         "Harmonic")

add_meta = function(result) {
  meta_result = lapply(result, function(x) {
    lapply(x, function(y) {
      y = y[c(1, 4, 5)]
      add_meta_result = c(3*min(y), 
                          exp(1)*prod(y)^(1/3),
                          exp(1)*log(3)*3/sum(1/y))
      return(c(y, add_meta_result))
    })
  })
  return(meta_result)
}

plot_meta = function(mean_power, legend_name, compare_methods, label_names) {
  mean_power = mean_power[c(1,4,2,5,3,6),]
  df_power = data.frame(mu_seq = rep(0:5, each = nrow(mean_power)),
                        power = as.vector(mean_power),
                        grp = rep(legend_name,
                                  ncol(mean_power)))
  p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
             aes(x = mu_seq, y = power, group = grp, fill = grp)) +
    geom_line(aes(linetype = grp, color = grp), size = 0.8) +
    geom_hline(yintercept=0.05) +
    geom_point(aes(shape = grp, color = grp), size = 2.5) +
    scale_shape_manual(values = c(21, 7, 22, 8, 3, 10), labels = label_names) +
    scale_color_manual(values = c("red1", "yellow4", "deepskyblue2", "tan1", "purple", "deeppink"),
                       labels = label_names) +
    scale_fill_manual(values = c("red1", "yellow4", "deepskyblue2", "tan1", "purple", "deeppink"),
                      labels = label_names) +
    scale_linetype_manual(values = rep(c("dashed", "solid"),3), labels = label_names) +
    theme(legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = c("none"), legend.text = element_text(size = 8)) +
    xlab("Scale of treatment effect") + ylab("power") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
  plot(p)
  return(p)
} 

mode = "dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_cauchy"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_dense_weak"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_sparse_strong"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_pos_strong"
load(paste("result/",mode,".Rdata", sep = ""))
meta_result = add_meta(result)
mean_power = sapply(meta_result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_meta(mean_power = mean_power, legend_name = legend_name,
              compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/meta_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)


########################### Figure 14 ###################################
mode = "linear_cauchy_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); 
legend_name[c(6, 8, 9)] = c("i-Wilcoxon-robust_old", "i-Wilcoxon-original",
                            "i-Wilcoxon-robust")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon-robust", "linear-CATE-test",
                    "i-Wilcoxon-original", "i-Wilcoxon-robust")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, 
            legend_name = legend_name, pos = c(0.25, 0.85))
ggsave(filename = paste("figure/interactive_",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)


#########################################################################
########################### small sample size  ##########################
#########################################################################


#########################################################################
########################### Figure 15 ###################################
#########################################################################
add_meta_small = function(result1, result2, robust = TRUE) {
  result = list()
  if(robust) {
    for(i in 1:length(result1)) {
      result[[i]] = list()
      for(j in 1:length(result1[[i]])){
        result[[i]][[j]] = c(result1[[i]][[j]][c("CovAdj-Wilcoxon-robust","linear-CATE-test","i-Wilcoxon-robust-signedA")],
                             3*min(result2[[i]][[j]][c(1,4,5)]),
                             4*min(result1[[i]][[j]]["i-Wilcoxon-robust-signedA"], c(result2[[i]][[j]][c(1,4,5)])))
      }
    }
  } else {
    for(i in 1:length(result1)) {
      result[[i]] = list()
      for(j in 1:length(result1[[i]])){
        result[[i]][[j]] = c(result1[[i]][[j]][c("CovAdj-Wilcoxon-quadratic","linear-CATE-test","i-Wilcoxon-quadratic-signedA")],
                             3*min(result2[[i]][[j]][c(1,4,5)]),
                             4*min(result1[[i]][[j]]["i-Wilcoxon-quadratic-signedA"], c(result2[[i]][[j]][c(1,4,5)])))
      }
    }
  }
  return(result)
}

mode = "linear_control_skewed"
load(paste("result/",mode,"_updated_small.Rdata", sep = ""))
i_result = result
load(paste("result/",mode,"_permutation.Rdata", sep = ""))
meta_result = result
result = add_meta_small(i_result, meta_result)
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
colnames(mean_power) = 0:5
legend_name = c("CovAdj-Wilcoxon", "linear-CATE-test",
                "i-Wilcoxon", "Wilcoxon-Bonferroni", "i-Wilcoxon-Bonferroni")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon", "Wilcoxon-Bonferroni", "i-Wilcoxon-Bonferroni")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none", legend_size = 6)
ggsave(filename = paste("figure/interactive_",mode,"_small.eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "linear_cauchy"
load(paste("result/",mode,"_updated_small.Rdata", sep = ""))
i_result = result
load(paste("result/",mode,"_permutation.Rdata", sep = ""))
meta_result = result
result = add_meta_small(i_result, meta_result)
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
colnames(mean_power) = 0:5
legend_name = c("CovAdj-Wilcoxon", "linear-CATE-test",
                "i-Wilcoxon", "Wilcoxon-Bonferroni", "i-Wilcoxon-Bonferroni")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon", "Wilcoxon-Bonferroni", "i-Wilcoxon-Bonferroni")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none", legend_size = 6)
ggsave(filename = paste("figure/interactive_",mode,"_small.eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "quadratic"
load(paste("result/",mode,"_updated_small.Rdata", sep = ""))
i_result = result
load(paste("result/",mode,"_permutation.Rdata", sep = ""))
meta_result = result
result = add_meta_small(i_result, meta_result, robust = FALSE)
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
colnames(mean_power) = 0:5
legend_name = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon", "Wilcoxon-Bonferroni", "i-Wilcoxon-Bonferroni")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("CovAdj-Wilcoxon", "linear-CATE-test",
                    "i-Wilcoxon", "Wilcoxon-Bonferroni", "i-Wilcoxon-Bonferroni")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none", legend_size = 6)
ggsave(filename = paste("figure/interactive_",mode,"_small.eps", sep = ""),
       plot = p, width = 4, height = 3.6)



########################### Figure 16 ###################################
mode = "linear_cauchy_oracle_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]]); legend_name[1] = "i-Wilcoxon-robust"
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon-robust", "i-Wilcoxon-oracle")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = c(0.25, 0.85), legend_size = 10)
ggsave(filename = "figure/oracle_power_cauchy.eps",
       plot = p, width = 4, height = 3.6)


#########################################################################
########################### Figure 17 ###################################
#########################################################################
plot_df_small = function(mean_power, compare_methods, label_names) {
  legend_name = factor(methods_var_Wilcoxon, levels = methods_var_Wilcoxon)
  df_power = data.frame(mu_seq = rep(as.numeric(colnames(mean_power)), each = nrow(mean_power)),
                        power = as.vector(mean_power),
                        grp = rep(legend_name,
                                  ncol(mean_power)))
  p = ggplot(data = subset(df_power, (grp %in% compare_methods)),
             aes(x = mu_seq, y = power, group = grp, fill = grp)) +
    geom_line(aes(linetype = grp, color = grp), size = 0.8) +
    geom_hline(yintercept=0.05) +
    geom_point(aes(shape = grp, color = grp), size = 2.5) +
    scale_shape_manual(values = c(21, 22, 3), labels = label_names) +
    scale_color_manual(values = c("red1", "deepskyblue2", "purple"), labels = label_names) +
    scale_fill_manual(values = c("red1", "deepskyblue2", "purple"), labels = label_names) +
    scale_linetype_manual(values = c("solid", "42", "44"), labels = label_names) +
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

compare_methods = c("R(X)",  "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
                    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
label_names = expression(R(X), 
                         paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|"),
                         S %.% (paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|")))

mode = "dense_weak_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_cauchy_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "dense_weak_control_skewed_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_smallN_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_smallN_cauchy_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "sparse_strong_smallN_control_skewed_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_dense_weak_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_sparse_strong_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_pos_strong_smallN_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df_small(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/",mode,".eps", sep = ""),
       plot = p, width = 4, height = 3.6)



############################ rebuttal Figure 1 ####################################
mode = "linear_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]])
legend_name[c(5,8)] = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,"_compare.eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "linear_control_skewed_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]])
legend_name[c(6,9)] = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,"_compare.eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "linear_cauchy_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]])
legend_name[c(6,9)] = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,"_compare.eps", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "quadratic_updated"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
legend_name = names(result[[1]][[1]])
legend_name[c(7,10)] = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
legend_name = factor(legend_name, levels = legend_name)
compare_methods = c("i-Wilcoxon (original)", "i-Wilcoxon (weighted)")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name, pos = "none")
ggsave(filename = paste("figure/interactive_",mode,"_compare.eps", sep = ""),
       plot = p, width = 4, height = 3.6)





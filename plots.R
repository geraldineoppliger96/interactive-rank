dir.create("figure")
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
          legend.position = c(0.8,0.52), legend.text = element_text(size = 6)) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    xlab("Scale of treatment effect") + ylab("power") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
  plot(p) 
  return(p)
}


############## Figure 2
mode = "linear_both"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
temp_name = names(result[[1]][[1]])
temp_name[c(2, 6, 5, 4)] = c("CovAdj-Wilcoxon", "interactive",
                             "interactive-robust", "interactive-quadratic")
legend_name = factor(temp_name, levels = temp_name[c(1,2,3,6,5,4)])
compare_methods = c("CovAdj-Wilcoxon", "CATE-test",
                    "interactive")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


mode = "blind_linear_both"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
temp_name = names(result[[1]][[1]])
temp_name[c(2, 5, 4)] = c("CovAdj-Wilcoxon", "interactive",
                             "interactive-robust")
legend_name = factor(temp_name, levels = temp_name[c(1,2,3,5,4)])
compare_methods = c("CovAdj-Wilcoxon", "CATE-test",
                    "interactive")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


############## Figure 3
mode = "linear_both_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                    byrow = TRUE, nrow = length(x)))})
temp_name = names(result[[1]][[1]])
temp_name[c(2, 6, 5, 4)] = c("CovAdj-Wilcoxon-robust", "interactive-original",
                             "interactive-robust", "interactive-quadratic")
legend_name = factor(temp_name, levels = temp_name[c(1,2,3,6,5,4)])
compare_methods = c("CovAdj-Wilcoxon-robust", "CATE-test",
                    "interactive-original", "interactive-robust")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


############## Figure 4
mode = "linear_pos"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(result, function(x){colMeans(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
sd_power = sapply(result, function(x){colSes(matrix(unlist(x) < alpha,
                                                        byrow = TRUE, nrow = length(x)))})
temp_name = names(result[[1]][[1]])
temp_name[c(2, 6, 5, 4)] = c("CovAdj-Wilcoxon-quadratic", "interactive-original",
                              "interactive-robust", "interactive-quadratic")
legend_name = factor(temp_name, levels = temp_name[c(1,2,3,6,5,4)])
compare_methods = c("CovAdj-Wilcoxon-quadratic", "CATE-test",
                    "interactive-robust", "interactive-quadratic")
p = plot_df(mean_power = mean_power, compare_methods = compare_methods,
            legend_name = legend_name)
ggsave(filename = paste("figure/interactive_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)



############### Figure 5(a)
mode = "varyx_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
mode = "cauchy_varyx_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power_cauchy = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R-Gaussian", "R-Cauchy", "R(X, 1-A)-Gaussian", "R(X, 1-A)-Cauchy")

legend_name = factor(compare_methods, levels = compare_methods)
df_combine = data.frame(mu_seq = rep(as.numeric(colnames(mean_power)),
                                  each = 4),
                     power = as.vector(rbind(mean_power[c(1,23),],
                                             mean_power_cauchy[c(1,23),])[c(1,3,2,4),]),
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


############### Figure 5(b)
mode = "varyx_big_more"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R", "R(X, 1-A)")

temp_name = methods_unpair
temp_name[c(1, 23, 20, 8, 15)] =
  c("R",
    "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
    "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
legend_name = factor(temp_name, levels = c(temp_name[c(1, 23, 20, 8, 15)],
                                           temp_name[-c(1, 23, 20, 8, 15)]))
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
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = c("red1", "green3")) +
  scale_fill_manual(values = c("red1", "green3")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
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


############### Figure 5(c)
mode = "varyx_small_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
compare_methods = c("R", "R(X, 1-A)")

temp_name = methods_unpair
temp_name[c(1, 23, 20, 8, 15)] =
  c("R",
    "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
    "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
legend_name = factor(temp_name, levels = c(temp_name[c(1, 23, 20, 8, 15)],
                                           temp_name[-c(1, 23, 20, 8, 15)]))
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
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = c("red1", "green3")) +
  scale_fill_manual(values = c("red1", "green3")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
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


############### Figure 6(a)
mode = "varyx_small_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R(X, 1-A)",  "R - hat(R)(X, 1 - A)")
  
temp_name = methods_unpair
temp_name[c(1, 23, 20, 8, 15)] =
  c("R",
    "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
    "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
legend_name = factor(temp_name, levels = c(temp_name[c(1, 23, 20, 8, 15)],
                                           temp_name[-c(1, 23, 20, 8, 15)]))
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


mode = "varyx_big_more_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R(X, 1-A)",  "R - hat(R)(X, 1 - A)")

temp_name = methods_unpair
temp_name[c(1, 23, 20, 8, 15)] =
  c("R",
    "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
    "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
legend_name = factor(temp_name, levels = c(temp_name[c(1, 23, 20, 8, 15)],
                                           temp_name[-c(1, 23, 20, 8, 15)]))
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


################ Figure 7
plot_df = function(mean_power, compare_methods, label_names) {
  temp_name = methods_unpair
  temp_name[c(1, 23, 20, 8, 15)] =
    c("R",
      "R(X, 1-A)",  "R - hat(R)(X, 1 - A)",
      "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
      "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
  legend_name = factor(temp_name, levels = c(temp_name[c(1, 23, 20, 8, 15)],
                                             temp_name[-c(1, 23, 20, 8, 15)]))
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

mode = "varyx_big_more"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "cauchy_varyx_big_more"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "varyx_big_more_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


################### Figure 8
mode = "varyx_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "cauchy_varyx_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "varyx_small_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "pos_big"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_big"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "both_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
names = expression(R(X, 1-A), R - hat(R)(X, 1 - A))
compare_methods = c("R",  "R - hat(R)(X, 1 - A)", "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|")
label_names = expression(R, R - hat(R)(X, 1 - A),
                         abs( R - hat(R)(X, 1 - A)) - abs( R - hat(R)(X, A)))
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/ranks_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)


################### Figure 9
compare_methods = c("R",  "R - hat(R)(X, 1 - A)",
                    "|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|",
                    "S(|R - hat(R)(X, 1 - A)| - |R - hat(R)(X, A)|)")
label_names = expression(R, R - hat(R)(X, 1 - A), 
                         paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|"),
                         S %.% (paste("|", R - hat(R)(X, 1 - A), "|") - paste("|", R - hat(R)(X, A), "|")))

mode = "both_small"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "varyx_big_more"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "varyx_big_more_control_skewed"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)

mode = "cauchy_varyx_big_more"
load(paste("result/",mode,".Rdata", sep = ""))
mean_power = sapply(rejection, function(x){
  colMeans(matrix(unlist(x) < alpha, byrow = TRUE, nrow = length(x)), na.rm = TRUE)})
p = plot_df(mean_power = mean_power, compare_methods = compare_methods, label_names = label_names)
ggsave(filename = paste("figure/signed_",mode,".png", sep = ""),
       plot = p, width = 4, height = 3.6)







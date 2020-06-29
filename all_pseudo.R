batch_pseudo_wilcoxon = function(dat, alg_type, n_permute, paired = FALSE) {
  W_obs = stat_generator(dat, alg_type = alg_type); print(W_obs)
  W_permute = foreach(i = 1:n_permute, .combine = cbind,
                      .export = c("stat_generator"),
                      .packages = c("foreach", "randomForest")) %dopar% {
                        dat_permute = dat; 
                        if (paired) {
                          temp = base::sample(dat$Z[1:(n/2)])
                          dat_permute$Z = c(temp, 1 - temp)
                        } else {
                          dat_permute$Z = base::sample(dat_permute$Z)
                        }
                       stat_generator(dat_permute, alg_type = alg_type)
                      }
  pval = rowMeans(apply(W_permute, 2, function(x) {x > W_obs}))
  return(pval)
}

if(0){
temp_name = c("covadj", "diff-in-error", "diff", "error")
W_dist = data.frame(grp = rep(factor(temp_name, levels = temp_name), each = 100),
                    w_dist = c(W_cov_permute, as.vector(t(W_permute[c(2,14,8),]))))
p = ggplot(W_dist, aes(x=grp, y=w_dist)) + # geom_boxplot(outlier.shape=NA) +
  geom_violin() +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=12)) +
  ylab("Null distribution") + ylim(-1100, 1350) +
  geom_segment(aes(x = 0.5, y = w_cov, xend = 1.5, yend = w_cov), size=1.2, color = "red") +
  geom_segment(aes(x = 1.5, y = W_obs[2], xend = 2.5, yend = W_obs[2]), size=1.2, color = "red") +
  geom_segment(aes(x = 2.5, y = W_obs[14], xend = 3.5, yend = W_obs[14]), size=1.2, color = "red") +
  geom_segment(aes(x = 3.5, y = W_obs[8], xend = 4.5, yend = W_obs[8]), size=1.2, color = "red") 
plot(p)
ggsave(filename = paste("figure/unpair_null_varyx_big_more_cauchy.png", sep = ""),
       plot = p, width = 4, height = 3.6)

ind = dat$X.3 > 1.5 #(1 - sin(dat$X.3)) > 0.7
y = E
x = (1 - abs(sin(dat$X.3))) > 0.5
# y = e
# x = (dat$Z == 1)
df = data.frame(y = y, x = x)
ggplot(df, aes(x=y, color=x, fill=x)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+xlim(-4,4) + 
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
  #            linetype="dashed")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(title="Weight histogram plot",x="Weight(kg)", y = "Density")+
  theme_classic()

}

stat_generator = function(dat, alg_type, max_nodes = NULL) {
  W = vector(length = 18)
  names(W) = c("rank-pseudo", "value-pseudo",
               "rank-signed-pseudo1", "value-signed-pseudo1",
               "rank-signed-pseudo2", "value-signed-pseudo2",
               "rank-signed-pseudo3", "value-signed-pseudo3",
               "rank-mask", "value-mask",
               "rank-false", "value-false", "mann-false",
               "rank-aver", "value-aver", "mann-aver",
               "rank-diff", "value-diff")
  
  if (alg_type == "linear") {
    model_base = lm(Y ~ ., data = dat)
    model_mask = lm(Y ~ .-Z, data = dat)
  } else if (alg_type == "RF") {
    model_base = randomForest(Y ~ ., data = dat, maxnodes = max_nodes)
    model_mask = randomForest(Y ~ .-Z, data = dat, maxnodes = max_nodes)
  }

  yhat_true = predict(model_base); yhat_mask = predict(model_mask)
  
  dat_false = dat; dat_false$Z = 1 - dat_false$Z
  yhat_false = predict(model_base, newdata = dat_false)
  
  E = abs(yhat_false - dat$Y) - abs(yhat_true - dat$Y)
  W["rank-pseudo"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-pseudo"] = sum(E)
  
  E = (abs(yhat_false - dat$Y) - abs(yhat_true - dat$Y))*
    sign((yhat_true - yhat_false)*(2*dat$Z - 1))
  W["rank-signed-pseudo1"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-signed-pseudo1"] = sum(E)
  
  E = (abs(yhat_false - dat$Y) - abs(yhat_true - dat$Y))*
    sign((dat$Y - yhat_false)*(2*dat$Z - 1))
  W["rank-signed-pseudo2"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-signed-pseudo2"] = sum(E)
  
  comb_ind = (dat$Y - yhat_false)*(2*dat$Z - 1) >= 0 |
    (yhat_true - yhat_false)*(2*dat$Z - 1) >= 0
  E = (abs(yhat_false - dat$Y) - abs(yhat_true - dat$Y))*
    (2*comb_ind - 1)
  W["rank-signed-pseudo3"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-signed-pseudo3"] = sum(E)
  
  E = abs(yhat_mask - dat$Y) - abs(yhat_true - dat$Y) 
  W["rank-mask"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-mask"] = sum(E)
  
  E = (dat$Y - yhat_false)*(2*dat$Z - 1)
  W["rank-false"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-false"] = sum(E)
  W["mann-false"] = sum(rank(dat$Y - yhat_false)*(2*dat$Z - 1))
  
  E = (2*dat$Y - yhat_false - yhat_true)*(2*dat$Z - 1)
  W["rank-aver"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-aver"] = sum(E)
  W["mann-aver"] = sum(rank(2*dat$Y - yhat_false - yhat_true)*(2*dat$Z - 1))
  
  E = (yhat_true - yhat_false)*(2*dat$Z - 1)
  #W["rank-diff"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["rank-diff"] = sum((2*(E > 0) - 1)*rank(abs(E)))
  W["value-diff"] = sum(E)
  return(W)
}



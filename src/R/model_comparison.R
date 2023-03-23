# Estimate the difference in the models expected predictive accuracy 
# using the difference in elpd_loo or elpd_waic
rm(list=ls())

library(loo)

threshold.seropositive <- 64  
n.age.classes <- 40
m0 <- 0.9

model_df <- data.frame(model_num = c("1","2"),
                       model_abb = c("mab_const","mab_rw"))
model_df$ll = rep(NA,times=nrow(model_df))

loo_l = vector(mode = "list", length = nrow(model_df))
waic_l = vector(mode = "list", length = nrow(model_df))

# fill in df
for(i in 1:nrow(model_df)){
  model_abb = model_df$model_abb[i]
  fit = readRDS(file=paste0("results/stan_fit/fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".rds"))
  
  log_lik = loo::extract_log_lik(fit, parameter_name = "log_likelihood", merge_chains = TRUE)
  
  model_df$ll[i] = sum(apply(log_lik,2,mean))
  
  loo_l[[i]] = loo(fit,pars="log_likelihood")
  # loo = loo(log_lik)
  
  waic_l[[i]] = waic(log_lik)
  
  pdf(paste0("figures/loo_plot_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".pdf"))
  plot(loo(log_lik),label_points = T)
  dev.off()
}

# Likelihood
file_name_ll = paste0("results/model_comparison/ll_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".csv")
write.csv(model_df,file=file_name_ll,row.names=FALSE)

# LOO
loo_compare(loo_l[[1]], loo_l[[2]])
loo_res = loo_compare(loo_l[[1]], loo_l[[2]])
file_name_loo = paste0("results/model_comparison/loo_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".csv")
write.csv(loo_res,file=file_name_loo,row.names=TRUE)

# WAIC
loo_compare(waic_l[[1]], waic_l[[2]])
waic_res = loo_compare(waic_l[[1]], waic_l[[2]])
file_name_waic = paste0("results/model_comparison/waic_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".csv")
write.csv(waic_res,file=file_name_waic,row.names=TRUE)


# Comparing best and 2nd best
loo_compare(loo_l[[2]],loo_l[[1]])
z_1 = abs(loo_compare(loo_l[[2]],loo_l[[1]])[2,1])/loo_compare(loo_l[[2]],loo_l[[1]])[2,2]
p_1 = pnorm(z_1,lower.tail = FALSE)
p_1



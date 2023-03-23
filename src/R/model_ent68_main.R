rm(list=ls())

library("rstan")
library("bayesplot")
library("RColorBrewer")
library("dplyr")
library("binom")

source("src/R/functions_model_ent68.R")

# Set model to run. Options are:
# "mab_const", "mab_rw"
model_abb = "mab_const"

#----------------------#
#-- Startup messages --#
#----------------------#

# If you are using rstan locally on a multicore machine and have plenty 
# of RAM to estimate your model in parallel, at this point execute
# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())

# To avoid recompilation of unchanged Stan programs, we recommend calling
# Allows you to automatically save a bare version of a compiled Stan program 
# to the hard disk so that it does not need to be recompiled (unless you 
# change it).
rstan_options(auto_write = TRUE)

#----------#
#-- Data --#
#----------#

threshold.seropositive <- 64
n.age.classes <- 40
m0 <- 0.9 # proportion seropositive at birth (1.0 for 1:16 and 0.9 for 1:64)

ent68_dat = get_data_for_stan_model_mab(threshold.seropositive,n.age.classes,m0)

#-------------------#
#-- Model fitting --#
#-------------------#

# Stan program
model_name_str = get_model_name_str(model_abb)

model <- stan_model(file=paste0("src/stan/",model_name_str))

# Call to fitting
fit <- sampling(model,                    # Stan program
                data = ent68_dat,         # named list of data
                chains = 4,               # number of Markov chains
                warmup = 2500,            # number of warmup iterations per chain (3000)
                iter = 10000,             # total number of iterations per chain (10000)
                refresh = 1000,           # show progress every 'refresh' iterations
                init = set_init_f(model_abb),
                control=list(adapt_delta=0.99,
                             max_treedepth=12)
)

# Save results in file
saveRDS(fit,file=paste0("results/stan_fit/fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".rds"))
# fit = readRDS(file=paste0("results/stan_fit/fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".rds"))


# Exploring the results

# library(shinystan)
# launch_shinystan(fit)

library(posterior)
sum_df <- summarise_draws(fit)
sum(sum_df$rhat[!is.na(sum_df$rhat)] > 1.01) == 0 
sum(sum_df$ess_tail[!is.na(sum_df$ess_tail)] < 400) == 0 
sum(sum_df$ess_bulk[!is.na(sum_df$ess_bulk)] < 400) == 0

# The object fit, returned from function stan is an S4 object of class stanfit
# fit: provides a summary for the parameter of the model as well as the log-posterior 
# with name lp__
print(fit)


#-----------------------#
#-- Analyse model fit --#
#-----------------------#

model_pars = get_model_parameters(fit)

### return an array of three dimensions: iterations, chains, parameters 
# permuted, a logical scalar indicating whether the draws after the warmup period in 
# each chain should be permuted and merged. If FALSE, the original order is kept.
a <- extract(fit, permuted = FALSE) 

median(a[,"chain:1","omega"])
median(a[,"chain:2","omega"])
median(a[,"chain:3","omega"])
median(a[,"chain:4","omega"])

#-----------------#
#-- Plot traces --#
#-----------------#

traces.colors <- rev(color_scheme_get("viridisA", c(2,3,4,5)))

op=par(mfrow=c(6,6),mar=c(1, 2, 2.5, 0.5) + 0.1,oma=c(2,1,1,0))
for(i in 1:length(model_pars)){
  param = model_pars[i]
  min.y = min(a[,,param])
  max.y = max(a[,,param])

  plot(a[,"chain:1",param],type="l",ylim=c(min.y,max.y),
       xlab="",ylab="",col=traces.colors[[1]]) # extract one chain, one parameter
  title(param, line = 0.5)
  lines(a[,"chain:2",param],col=traces.colors[[2]])
  lines(a[,"chain:3",param],col=traces.colors[[3]])
  lines(a[,"chain:4",param],col=traces.colors[[4]])
}
par(op)


#--------------------#
#-- Plot densities --#
#--------------------#

op=par(mfrow=c(6,6),mar=c(1, 2, 2.5, 0.5) + 0.1,oma=c(2,1,1,0))
for(i in 1:length(model_pars)){
  param = model_pars[i]

  plot(density(a[,"chain:1",param]),col=traces.colors[[1]],xlab="",ylab="",main="")
  title(param, line = 0.5)
  lines(density(a[,"chain:2",param]),col=traces.colors[[2]])
  lines(density(a[,"chain:3",param]),col=traces.colors[[3]])
  lines(density(a[,"chain:4",param]),col=traces.colors[[4]])
}
par(op)


#----------------------------#
#-- Plot model fit to data --#
#----------------------------#

# Extract fit summary (mean, se, sd...)
fit.summary <- data.frame(summary(fit)$summary)

fit_cols = brewer.pal(n=8,"Dark2")[1:2]

years_v = c(ent68_dat$year1,ent68_dat$year2,ent68_dat$year3)

pr_pos_1_2006.median <- c(m0,fit.summary$X50.[grep("pr_pos_1_2006",rownames(fit.summary))])
pr_pos_1_2006.low <- c(m0,fit.summary$X2.5.[grep("pr_pos_1_2006",rownames(fit.summary))])
pr_pos_1_2006.up <- c(m0,fit.summary$X97.5.[grep("pr_pos_1_2006",rownames(fit.summary))])

pr_pos_2006.median <- fit.summary$X50.[grep("pr_pos_2006",rownames(fit.summary))]
pr_pos_2006.low <- fit.summary$X2.5.[grep("pr_pos_2006",rownames(fit.summary))]
pr_pos_2006.up <- fit.summary$X97.5.[grep("pr_pos_2006",rownames(fit.summary))]

pr_pos_1_2011.median <- c(m0,fit.summary$X50.[grep("pr_pos_1_2011",rownames(fit.summary))])
pr_pos_1_2011.low <- c(m0,fit.summary$X2.5.[grep("pr_pos_1_2011",rownames(fit.summary))])
pr_pos_1_2011.up <- c(m0,fit.summary$X97.5.[grep("pr_pos_1_2011",rownames(fit.summary))])

pr_pos_2011.median <- fit.summary$X50.[grep("pr_pos_2011",rownames(fit.summary))]
pr_pos_2011.low <- fit.summary$X2.5.[grep("pr_pos_2011",rownames(fit.summary))]
pr_pos_2011.up <- fit.summary$X97.5.[grep("pr_pos_2011",rownames(fit.summary))]

pr_pos_1_2017.median <- c(m0,fit.summary$X50.[grep("pr_pos_1_2017",rownames(fit.summary))])
pr_pos_1_2017.low <- c(m0,fit.summary$X2.5.[grep("pr_pos_1_2017",rownames(fit.summary))])
pr_pos_1_2017.up <- c(m0,fit.summary$X97.5.[grep("pr_pos_1_2017",rownames(fit.summary))])

pr_pos_2017.median <- fit.summary$X50.[grep("pr_pos_2017",rownames(fit.summary))]
pr_pos_2017.low <- fit.summary$X2.5.[grep("pr_pos_2017",rownames(fit.summary))]
pr_pos_2017.up <- fit.summary$X97.5.[grep("pr_pos_2017",rownames(fit.summary))]

# DF
post_est_df = expand.grid(Decimal_age=c(0,ent68_dat$Decimal_age_under1,ent68_dat$Decimal_age),year=years_v)
# Add in data
post_est_df$n = c(c(0,ent68_dat$n1_2006,ent68_dat$n_2006),
                  c(0,ent68_dat$n1_2011,ent68_dat$n_2011),
                  c(0,ent68_dat$n1_2017,ent68_dat$n_2017))
post_est_df$z = c(c(0,ent68_dat$z1_2006,ent68_dat$z_2006),
                  c(0,ent68_dat$z1_2011,ent68_dat$z_2011),
                  c(0,ent68_dat$z1_2017,ent68_dat$z_2017))
# Add in 95% CI for the data
post_est_df = post_est_df %>%
  mutate(binom.confint(x=z,n=n,conf.level=0.95,methods="exact")[5:6] %>%
           rename(seroprev_lower=lower,seroprev_upper=upper))

# Add in posterior estimates
post_est_df$post_pr_pos.median = c(pr_pos_1_2006.median,pr_pos_2006.median,
                                   pr_pos_1_2011.median,pr_pos_2011.median,
                                   pr_pos_1_2017.median,pr_pos_2017.median)
post_est_df$post_pr_pos.low = c(pr_pos_1_2006.low,pr_pos_2006.low,
                                pr_pos_1_2011.low,pr_pos_2011.low,
                                pr_pos_1_2017.low,pr_pos_2017.low)
post_est_df$post_pr_pos.up = c(pr_pos_1_2006.up,pr_pos_2006.up,
                               pr_pos_1_2011.up,pr_pos_2011.up,
                               pr_pos_1_2017.up,pr_pos_2017.up)

# Plot
p <- ggplot(data=post_est_df,aes(x=Decimal_age)) +
  geom_pointrange(aes(y=z/n,ymin=seroprev_lower,ymax=seroprev_upper),
             col=fit_cols[1],fatten=1.5,shape=15) +
  geom_line(aes(y=post_pr_pos.median),col=fit_cols[2]) +
  geom_ribbon(aes(ymin=post_pr_pos.low,ymax=post_pr_pos.up),alpha=0.5,
              fill=fit_cols[2],color="transparent") +
  labs(title = "Model fit - Posterior estimates of seroprevalence",
       x = "Age (years)", y = "Seroprevalence") +
  ylim(0,1) +
  facet_wrap("year") +
  theme_bw() 

p = p + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                           breaks=c(0,1,2,5,10,20,40),
                           labels=c(0,1,2,5,10,20,40)) 
p

ggsave(filename = paste0("figures/model_fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_with_under1.pdf"),
       p,height=90,width=183,units="mm")

#--------------------#
#-- Log-likelihood --#
#--------------------#

ll <- fit.summary[grep("log_likelihood",row.names(fit.summary)),]
sum(ll$mean)

loo(fit,pars="log_likelihood")



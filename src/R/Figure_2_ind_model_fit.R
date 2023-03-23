rm(list=ls())

library(rstan)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(binom)

source("src/R/functions_model_ent68.R")
source("src/R/my_colors.R")

savefigure = TRUE

model_abb = "mab_const"
threshold.seropositive <- 64
n.age.classes <- 40
m0 <- 0.9

# Read stan object
fit = readRDS(file=paste0("results/stan_fit/fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".rds"))

# Summary fit
fit.summary <- data.frame(summary(fit)$summary)

# Extract posterior estimates 
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


# Read in data
ent68_dat = get_data_for_stan_model_mab(threshold.seropositive,n.age.classes,m0)
years_v = c(ent68_dat$year1,ent68_dat$year2,ent68_dat$year3)


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

if(threshold.seropositive==16){ model_col = cutoff_cols[1]}
if(threshold.seropositive==64){ model_col = cutoff_cols[2]}

p <- ggplot(data=post_est_df,aes(x=Decimal_age)) +
  geom_pointrange(aes(y=z/n,ymin=seroprev_lower,ymax=seroprev_upper),
                  col="gray",fatten=1.5,shape=15) +
  geom_line(aes(y=post_pr_pos.median),col=model_col) +
  geom_ribbon(aes(ymin=post_pr_pos.low,ymax=post_pr_pos.up),alpha=0.5,
              fill=model_col,color="transparent") +
  labs(#title = "Model fit - Posterior estimates of seroprevalence",
       x = "Age (years)", y = "Seroprevalence") +
  ylim(0,1) +
  facet_wrap("year") +
  theme_bw() 

p = p + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                           breaks=c(0,1,2,5,10,20,40),
                           labels=c(0,1,2,5,10,20,40)) 
p

saveRDS(p,file=paste0("figures/model_fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".rds"))

if(savefigure){
  ggsave(filename = paste0("figures/model_fit_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".pdf"),
         height=90,width=183,units="mm")
}



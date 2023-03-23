rm(list=ls())

library(rstan)
library(RColorBrewer)

source("src/R/functions_model_ent68.R")
source("src/R/my_colors.R")

savefigure = TRUE

model_abb = "mab_rw"
n.age.classes <- 40
m0_16 <- 1.0
m0_64 <- 0.9

# Read stan objects
fit_16 = readRDS(file=paste0("results/stan_fit/fit_",model_abb,"_threshold_",16,"_age_classes_",n.age.classes,"_m0_",m0_16,".rds"))
fit_64 = readRDS(file=paste0("results/stan_fit/fit_",model_abb,"_threshold_",64,"_age_classes_",n.age.classes,"_m0_",m0_64,".rds"))

# Summary fit
fit_16.summary <- data.frame(summary(fit_16)$summary)
fit_64.summary <- data.frame(summary(fit_64)$summary)

# Get model pars
model_pars = get_model_parameters(fit_16)

# Extract posterior estimates for proportion with mab
ma_1_2006.median_16 <- fit_16.summary$X50.[grep("ma_1_2006",rownames(fit_16.summary))]
ma_1_2006.low_16 <- fit_16.summary$X2.5.[grep("ma_1_2006",rownames(fit_16.summary))]
ma_1_2006.up_16 <- fit_16.summary$X97.5.[grep("ma_1_2006",rownames(fit_16.summary))]

ma_2006.median_16 <- fit_16.summary$X50.[grep("ma_2006",rownames(fit_16.summary))]
ma_2006.low_16 <- fit_16.summary$X2.5.[grep("ma_2006",rownames(fit_16.summary))]
ma_2006.up_16 <- fit_16.summary$X97.5.[grep("ma_2006",rownames(fit_16.summary))]

ma_1_2006.median_64 <- fit_64.summary$X50.[grep("ma_1_2006",rownames(fit_64.summary))]
ma_1_2006.low_64 <- fit_64.summary$X2.5.[grep("ma_1_2006",rownames(fit_64.summary))]
ma_1_2006.up_64 <- fit_64.summary$X97.5.[grep("ma_1_2006",rownames(fit_64.summary))]

ma_2006.median_64 <- fit_64.summary$X50.[grep("ma_2006",rownames(fit_64.summary))]
ma_2006.low_64 <- fit_64.summary$X2.5.[grep("ma_2006",rownames(fit_64.summary))]
ma_2006.up_64 <- fit_64.summary$X97.5.[grep("ma_2006",rownames(fit_64.summary))]


# Read in data
ent68_dat_16 = get_data_for_stan_model_mab(16,n.age.classes,m0_16)
ent68_dat_64 = get_data_for_stan_model_mab(64,n.age.classes,m0_64)
years_v = c(ent68_dat_16$year1,ent68_dat_16$year2,ent68_dat_16$year3)

# DF
post_est_df = data.frame(Decimal_age=c(0,ent68_dat_16$Decimal_age_under1,ent68_dat_16$Decimal_age))

# Add in posterior estimates
post_est_df$post_ma.median_16 = c(ent68_dat_16$m0,ma_1_2006.median_16,ma_2006.median_16)
post_est_df$post_ma.low_16 = c(ent68_dat_16$m0,ma_1_2006.low_16,ma_2006.low_16)
post_est_df$post_ma.up_16 = c(ent68_dat_16$m0,ma_1_2006.up_16,ma_2006.up_16)
post_est_df$post_ma.median_64 = c(ent68_dat_64$m0,ma_1_2006.median_64,ma_2006.median_64)
post_est_df$post_ma.low_64 = c(ent68_dat_64$m0,ma_1_2006.low_64,ma_2006.low_64)
post_est_df$post_ma.up_64 = c(ent68_dat_64$m0,ma_1_2006.up_64,ma_2006.up_64)

# Plot

# Remove ages > 20 yo
post_est_df = post_est_df[post_est_df$Decimal_age < 20,]

p <- ggplot(data=post_est_df,aes(x=Decimal_age)) +
  geom_line(aes(y=post_ma.median_16,color=cutoff_cols[1])) +
  geom_ribbon(aes(ymin=post_ma.low_16,ymax=post_ma.up_16),alpha=0.5,
              fill=cutoff_cols[1],color="transparent") +
  geom_line(aes(y=post_ma.median_64,color=cutoff_cols[2])) +
  geom_ribbon(aes(ymin=post_ma.low_64,ymax=post_ma.up_64),alpha=0.5,
              fill=cutoff_cols[2],color="transparent") +
  labs(x = "Age (years)", y = "Proportion with maternal antibodies") +
  ylim(0,1) +
  scale_color_manual(name="Cut-off",values = c("1:16"=cutoff_cols[1],"1:64"=cutoff_cols[2])) +
  theme_classic() 
  
p = p + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                           breaks=c(0,1,2,5,10,20),
                           labels=c(0,1,2,5,10,20)) +
  theme(legend.position = c(.75, .75)) 
p

saveRDS(p,file=paste0("figures/estimated_prop_w_mab_",model_abb,".rds"))

if(savefigure){
  ggsave(filename = paste0("figures/Figure_3_B_estimated_prop_w_mab.pdf"),
         height=90,width=90,units="mm")
}





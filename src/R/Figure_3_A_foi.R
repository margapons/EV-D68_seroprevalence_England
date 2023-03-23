rm(list=ls())

library(rstan)
library(RColorBrewer)

source("src/R/functions_model_ent68.R")
source("src/R/my_colors.R")

savefigure = FALSE

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

# Extract FOI estimates
foi_16_median = rev(fit_16.summary$X50.[grep("lambda",model_pars)])
foi_16_low = rev(fit_16.summary$X2.5.[grep("lambda",model_pars)])
foi_16_up = rev(fit_16.summary$X97.5.[grep("lambda",model_pars)])
foi_64_median = rev(fit_64.summary$X50.[grep("lambda",model_pars)])
foi_64_low = rev(fit_64.summary$X2.5.[grep("lambda",model_pars)])
foi_64_up = rev(fit_64.summary$X97.5.[grep("lambda",model_pars)])

# Years
ent68_dat = get_data_for_stan_model_mab(16,n.age.classes,m0_16)
min.year = ent68_dat$year1 - ent68_dat$N 
max.year = ent68_dat$year3 
x_axis = seq(min.year,max.year)

# Construct df for plot
foi_df <- data.frame(x_axis = x_axis,
                     foi_16_median = foi_16_median,
                     foi_16_low = foi_16_low,
                     foi_16_up = foi_16_up,
                     foi_64_median = foi_64_median,
                     foi_64_low = foi_64_low,
                     foi_64_up = foi_64_up
)
  

# Plot

survey_col = "gray" # color for arrows indicating sero-surveys

p <- ggplot(data=foi_df,aes(x=x_axis)) +
  geom_line(aes(y=foi_16_median,colour=cutoff_cols[1])) +
  geom_ribbon(aes(ymin=foi_16_low,ymax=foi_16_up),alpha=0.5,
              fill=cutoff_cols[1],color="transparent") +
  geom_line(aes(y=foi_64_median,colour=cutoff_cols[2])) +
  geom_ribbon(aes(ymin=foi_64_low,ymax=foi_64_up),alpha=0.5,
              fill=cutoff_cols[2],color="transparent") +
  geom_segment(aes(x = 2017, y = 1.25,
                   xend = 2017, yend = 1.1),
               arrow = arrow(length = unit(0.4, "cm")),col=survey_col) +
  geom_segment(aes(x = 2011, y = 1.25,
                   xend = 2011, yend = 1.1),
               arrow = arrow(length = unit(0.4, "cm")),col=survey_col) +
  geom_segment(aes(x = 2006, y = 1.25,
                   xend = 2006, yend = 1.1),
               arrow = arrow(length = unit(0.4, "cm")),col=survey_col) +
  labs(x = "Year", y = "Force of infection") +
  ylim(0,1.25) +
  scale_x_continuous(breaks=seq(1965,2020,by=5),
                     limits=c(1990,2017)) +
  scale_color_identity(name = "Cut-off",
                       breaks = c(cutoff_cols[1], cutoff_cols[2]),
                       labels = c("1:16", "1:64"),
                       guide = "legend") +
  theme_classic()

p = p + theme(legend.position = c(.25, .75)) 

p
    
saveRDS(p,file=paste0("figures/foi_",model_abb,".rds"))

if(savefigure){
  ggsave(filename = "figures/Figure_3_A_foi_mab_rw.pdf",p,height=90,width=90,units="mm")
}
  






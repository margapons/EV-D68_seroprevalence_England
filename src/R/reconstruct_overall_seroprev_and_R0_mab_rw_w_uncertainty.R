# Reconstruct seroprevalence and R0
# for the RW model using demographic data from England
# Note this implementation only works for the RW model
rm(list=ls())

library(patchwork)
library(dplyr)
library(ggplot2)

source("src/R/functions_model_ent68.R")
source("src/R/my_colors.R")

savefigure = TRUE

# Set seropositivity cut-off and number of age classes
model_abb = "mab_rw"
threshold.seropositive = 64
n.age.classes = 40  
m0 = 0.9

# England population size per age and year
age_year_lf_interp <- read.csv(file=paste0("data/population_by_age_England_0_to_",n.age.classes,"_interpolated.csv"),header=TRUE)
pop_size <- read.csv(file=paste0("data/population_by_age_England_0_to_",n.age.classes,"_interpolated.csv"),header=TRUE)

# Get estimated model parameters
fit = readRDS(file=paste0("results/stan_fit/fit_mab_rw",
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,".rds"))

# Get data
ent68_dat = get_data_for_stan_model_mab(threshold.seropositive,n.age.classes,m0)

# Set min and max year 
min.year = ent68_dat$year1 - ent68_dat$N 
max.year = ent68_dat$year3 

# DF for results
res_df <- expand.grid(Year=2006:2017,
                      Decimal_age=c(0,ent68_dat$Decimal_age_under1,ent68_dat$Decimal_age))

pop_size_years = unique(pop_size$Year)
res_df$Pop = rep(NA,times=nrow(res_df))

pop_size_function <- function(year,decimal_age){
  if(decimal_age == 0)
    return(0)
  else
    if(decimal_age > 0 & decimal_age < 1)
      return(pop_size$Pop[pop_size$Year==year & pop_size$Age == 0]/12)
    else
      return(pop_size$Pop[pop_size$Year==year & pop_size$Age == floor(decimal_age)])
}

for(i in 1:nrow(res_df)){
  res_df$Pop[i] = pop_size_function(res_df$Year[i],res_df$Decimal_age[i])
}

# Recover parameter estimates 
model_pars = get_model_parameters(fit)

# Sample size
sample_size = 500
total_sample_size = dim(fit)[1]*dim(fit)[2]
sample_idx_v = sample(x=seq(1,total_sample_size),size=sample_size,replace=F)

# Extract sample for parameter  omega
est_omega <- data.frame(sample_idx = sample_idx_v)
est_omega$Value = rstan::extract(fit, "omega")[[1]][sample_idx_v]

# Extract sample for parameters lambda's
lambda_v <- model_pars[grep("lambda",model_pars)]
for(i in 1:length(lambda_v)){
  # Extract sample for a single lambda
  est_lambda_single <- data.frame(sample_idx = sample_idx_v)
  est_lambda_single$Value = rstan::extract(fit, lambda_v[i])[[1]][sample_idx_v]
  est_lambda_single$Year = rev(seq(min.year,max.year))[i]
  
  # rbind
  if(i == 1){
    est_lambda = est_lambda_single
  }else{
    est_lambda = rbind(est_lambda,est_lambda_single)    
  }
}

# DF for results
sample_id_df = as.data.frame(sample_idx_v)
names(sample_id_df) = "sample_id"
res_df = merge(res_df,sample_id_df)

# Add in new vars
res_df$prob_seroneg = rep(NA,times=nrow(res_df))

# Loop over samples
for(sample_idx in 1:nrow(sample_id_df)){
  sample = sample_id_df$sample_id[sample_idx]
  
  # Parameter estimates for that sample
  omega = est_omega$Value[est_omega$sample_idx == sample]
  lambda = est_lambda$Value[est_lambda$sample_idx == sample]
  
  # Loop over years of interest
  for(year in 2006:2017){
    pr_neg = rep(NA,times=ent68_dat$N)
    pr_neg_1 = rep(NA,times=ent68_dat$N1)
    lambda_sum = rep(NA,times=ent68_dat$N)
    lambda_sum_1 = rep(NA,times=ent68_dat$N1)
    year3 = ent68_dat$year3
    year2 = year
    m0 = ent68_dat$m0
    
    # Loop over age classes >1
    for(i in 1:ent68_dat$N){
      pr_neg[i] = (1-m0)*exp(-sum(lambda[(year3-year2+1):(year3-year2+1+i)])-0.5*lambda[year3-year2+1])
      # Terms j<i
      for(j in 0:(i-1)){
        lambda_sum[j+1] = 0
        if(j<(i-1)){
          for(k in (j+1):(i-1)){
            lambda_sum[j+1] = lambda_sum[j+1] + lambda[year3-year2+i+1-k]
          }
        }
        lambda_sum[j+1] = lambda_sum[j+1] + 0.5*lambda[year3-year2+i+1-i]
        pr_neg[i] = pr_neg[i] + m0*exp(-j*omega)*(omega*(exp(-omega)-exp(-lambda[year3-year2+i+1-j]))/(lambda[year3-year2+i+1-j]-omega))*exp(-lambda_sum[j+1])
      }
      # Term j=i
      pr_neg[i] = pr_neg[i] + m0*exp(-i*omega)*(omega*(exp(-0.5*omega)-exp(-0.5*lambda[year3-year2+i+1-i]))/(lambda[year3-year2+i+1-i]-omega)) 
    }

    # Loop over age classes <1
    for(i in 1:ent68_dat$N1){
      pr_neg_1[i] = (1-m0)*exp(-(i-0.5)/12.0*lambda[year3-year2+1]);
      # Terms j<i
      if(i>1){
        for(j in 1:(i-1)){
          lambda_sum_1[j] = 0
          if(j<(i-1)){
            for(k in (j+1):(i-1)){
              lambda_sum_1[j] = lambda_sum_1[j] + 1.0/12*lambda[year3-year2+i+1-i]
            }
          }
          lambda_sum_1[j] = lambda_sum_1[j] + 0.5/12*lambda[year3-year2+i+1-i]
          pr_neg_1[i] = pr_neg_1[i] + m0*exp(-(j-1)/12.0*omega)*(omega*(exp(-1.0/12*omega)-exp(-1.0/12*lambda[year3-year2+i+1-i]))/(lambda[year3-year2+i+1-i]-omega))*exp(-lambda_sum_1[j])
        }
      }
      # Term j=i
      pr_neg_1[i] = pr_neg_1[i] + m0*exp(-(i-1)/12.0*omega)*(omega*(exp(-0.5/12*omega)-exp(-0.5/12*lambda[year3-year2+i+1-i]))/(lambda[year3-year2+i+1-i]-omega))
    }
    
    # Copy results for a sample and year into res_df
    res_df$prob_seroneg[res_df$sample_id == sample & res_df$Year == year & res_df$Decimal_age %in% ent68_dat$Decimal_age] = pr_neg
    res_df$prob_seroneg[res_df$sample_id == sample & res_df$Year == year & res_df$Decimal_age %in% ent68_dat$Decimal_age_under1] = pr_neg_1
    res_df$prob_seroneg[res_df$sample_id == sample & res_df$Year == year & res_df$Decimal_age == 0] = 1 - ent68_dat$m0
  }
  
}

# Add in var
res_df$prob_seropos = 1 - res_df$prob_seroneg


# Seroprevalence by age and year - compute median and 95% CrI  
seroprevalence_df = res_df %>%
  group_by(Year,Decimal_age) %>%
  summarise(seroprev_median = median(x=prob_seropos,na.rm=T),
            seroprev_low = quantile(x=prob_seropos,probs=.025,na.rm=T),
            seroprev_up = quantile(x=prob_seropos,probs=.975,na.rm=T)
            )

# Plot
p1 <- ggplot(data=seroprevalence_df[seroprevalence_df$Year %in% seq(2006,2017),],
             aes(x=Decimal_age, y = seroprev_median, color=factor(Year))) +
  geom_line() +
  geom_ribbon(aes(ymin=seroprev_low,ymax=seroprev_up,fill=factor(Year)),alpha=0.3,color="transparent") +
  ylim(0,1) +
  xlim(0,n.age.classes) +
  labs(x = "Age (years)", y = "Estimated seroprevalence",
       color="Year", fill="Year") +
  scale_fill_manual(values=prog_2006_2017) +
  scale_color_manual(values=prog_2006_2017) +
  theme_classic()

saveRDS(p1,file=paste0("figures/seroprevalence_by_age_and_year_",model_abb,
                      "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".rds"))

if(savefigure){
  ggsave(filename = paste0("figures/seroprevalence_by_age_and_year_",model_abb,
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".pdf"),
         p1,height=90,width=100,units="mm")
}

# Age-weighted seroprevalence for each year (w=pop, x=prob_seropos) 
overall_res = res_df %>%
  group_by(Year,sample_id) %>%
  summarise(n = n(),
            overall_seroprev = weighted.mean(x=prob_seropos,w=Pop)) %>%
  group_by(Year) %>%
  summarise(n = n(),
            overall_seroprev_median = median(x=overall_seroprev,na.rm=T),
            overall_seroprev_low = quantile(x=overall_seroprev,probs=.025,na.rm=T),
            overall_seroprev_up = quantile(x=overall_seroprev,probs=.975,na.rm=T))

# Adding in R0 using classical result R0 = 1 / (1 - Seroprev) from Anderson & May, not accounting for age-dependency
overall_res$R0_AM_median = 1 / (1 - overall_res$overall_seroprev_median)
overall_res$R0_AM_low = 1 / (1 - overall_res$overall_seroprev_low)
overall_res$R0_AM_up = 1 / (1 - overall_res$overall_seroprev_up)

# Save results in file
saveRDS(overall_res,file=paste0("results/overall_seroprevalence_and_R0_",model_abb,
                                "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".rds"))

# Plots
p2 <- ggplot(data=overall_res[overall_res$Year %in% seq(2006,2017),],
             aes(x=Year, y=overall_seroprev_median)) +
  geom_line() +
  geom_ribbon(aes(ymin=overall_seroprev_low,ymax=overall_seroprev_up),alpha=0.3,color="transparent") +
  ylim(0,1) +
  scale_x_continuous(breaks=seq(2006,2017,by=2)) +  
  labs(y = "Overall seroprevalence") +
  theme_bw()

p3 <- ggplot(data=overall_res[overall_res$Year %in% seq(2006,2017),],
             aes(x=Year, y=R0_AM_median)) +
  geom_line() +
  geom_ribbon(aes(ymin=R0_AM_low,ymax=R0_AM_up),alpha=0.3,color="transparent") +
  scale_x_continuous(breaks=seq(2006,2017,by=2)) +  
  labs(y = "R0") +
  theme_bw()

p4 <- p2 + p3
p4

if(savefigure){
  ggsave(filename = paste0("figures/overall_seroprevalence_and_R0_",model_abb,
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".pdf"),
         p4,height=90,width=183,units="mm")
}


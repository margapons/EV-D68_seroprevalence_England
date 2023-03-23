# Reconstruct number of new infections per year and age class 
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

# Add in vars 
# probability to be seronegative until previous year
res_df$prob_seroneg_year_minus_one = rep(NA,times=nrow(res_df))
# probability to seroconvert current year
res_df$prob_seroconvert_current_year = rep(NA,times=nrow(res_df))

for(i in 1:nrow(res_df)){
  year_i = res_df$Year[i]
  decimal_age_i = res_df$Decimal_age[i]
  sample_i = res_df$sample_id[i]
  
  # Parameter estimates
  omega = est_omega$Value[est_omega$sample_idx == sample_i]
  lambda = est_lambda$Value[est_lambda$sample_idx == sample_i]
  
  # Fill in probability to be seronegative
  if(year_i > 2006){ 
    if(decimal_age_i > 2){
      res_df$prob_seroneg_year_minus_one[i] = res_df$prob_seroneg[res_df$Year==year_i-1 & res_df$Decimal_age==decimal_age_i-1 & res_df$sample_id==sample_i]  
    }
    if(decimal_age_i == 1.5){
      res_df$prob_seroneg_year_minus_one[i] = ent68_dat$m0*omega*(exp(-0.5*omega)-exp(-0.5*lambda[ent68_dat$year3-year_i+2]))/(lambda[ent68_dat$year3-year_i+2]-omega) + (1 - ent68_dat$m0)*exp(-0.5*lambda[ent68_dat$year3-year_i+2])
    }
  }
  
  # Fill in probability to seroconvert
  res_df$prob_seroconvert_current_year[i] = 1 - exp(-0.5*(lambda[ent68_dat$year3-year_i+1] + lambda[ent68_dat$year3-year_i+2]))
}


# Add in var
res_df$new_infections = res_df$prob_seroneg_year_minus_one*res_df$prob_seroconvert_current_year*res_df$Pop

# Number of infections by age and year - compute median and 95% CrI of 
new_infections_df = res_df %>%
  group_by(Year,Decimal_age) %>%
  summarise(inf_median = median(x=new_infections,na.rm=T),
            inf_low = quantile(x=new_infections,probs=.025,na.rm=T),
            inf_up = quantile(x=new_infections,probs=.975,na.rm=T)
            )

# Save results in file
saveRDS(res_df,file=paste0("results/number_infections/number_infections_",model_abb,
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".rds"))

# Plot
p1 <- ggplot(data=new_infections_df[new_infections_df$Year %in% seq(2007,2017),],
             aes(x=Decimal_age, y = inf_median, color=factor(Year))) +
  geom_line() +
  # geom_ribbon(aes(ymin=inf_low,ymax=inf_up,fill=factor(Year)),alpha=0.3,color="transparent") +
  labs(x = "Age (years)", y = "Number of infections",
       color="Year", fill="Year") +
  scale_x_continuous(breaks=seq(1,11,by=1),limits=c(1,11)) +
  scale_color_manual(values=prog_2006_2017) +
  theme_classic()

saveRDS(p1,file=paste0("figures/number_infections_by_age_and_year_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",
                       m0,"_sample_size_",sample_size,".rds"))

if(savefigure){
  ggsave(filename = paste0("figures/number_infections_by_age_and_year_",model_abb,
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_sample_size_",sample_size,"_zoom.pdf"),
         p1,height=90,width=100,units="mm")
}

p2 <- ggplot(data=new_infections_df[new_infections_df$Year %in% seq(2007,2017),],
             aes(x=Decimal_age, y = inf_median, color=factor(Year))) +
  geom_line() +
  # geom_ribbon(aes(ymin=inf_low,ymax=inf_up,fill=factor(Year)),alpha=0.3,color="transparent") +
  labs(x = "Age (years)", y = "Number of infections",
       color="Year", fill="Year") +
  scale_x_continuous(breaks=seq(0,41,by=5),limits=c(1,41)) +
  scale_color_manual(values=prog_2006_2017) +
  theme_classic()

# saveRDS(p2,file=paste0("figures/number_infections_by_age_and_year_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,".rds"))

if(savefigure){
  ggsave(filename = paste0("figures/number_infections_by_age_and_year_",model_abb,
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".pdf"),
         p2,height=90,width=100,units="mm")
}

# Compute mean age at infection each year
res_df_v2 = res_df[!(res_df$Decimal_age < 1.5),]
mean_age_infection_df = res_df_v2 %>%
  group_by(Year,sample_id) %>%
  summarise(mean_age = weighted.mean(x=Decimal_age,w=new_infections))

mean_age_infection_df_v2 = mean_age_infection_df %>%
  group_by(Year) %>%
  summarise(median = median(x=mean_age,na.rm=T),
            low = quantile(x=mean_age,probs=.025,na.rm=T),
            up = quantile(x=mean_age,probs=.975,na.rm=T)
  )

saveRDS(mean_age_infection_df_v2,file=paste0("results/mean_age_at_infection_",
                                             model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,
                                             "_m0_",m0,"_sample_size_",sample_size,".rds"))

p3 <- ggplot(data=mean_age_infection_df_v2, aes(x=Year)) + 
  geom_pointrange(aes(y=median, ymin=low, ymax=up),
                  fatten=1.5,shape=10) +
  labs(x = "Year", y = "Mean age at infection") +
  scale_x_continuous(breaks=seq(2007,2017,by=2)) +
  theme_classic()

saveRDS(p3,file=paste0("figures/mean_age_at_infection_",model_abb,"_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,
                       "_m0_",m0,"_sample_size_",sample_size,".rds"))

if(savefigure){
  ggsave(filename = paste0("figures/mean_age_at_infection_",model_abb,
                           "_threshold_",threshold.seropositive,"_age_classes_",n.age.classes,"_m0_",m0,"_sample_size_",sample_size,".pdf"),
         p3,height=90,width=100,units="mm")
}



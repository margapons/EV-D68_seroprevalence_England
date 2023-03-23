rm(list=ls())

library(cowplot)

source("src/R/my_colors.R")

savefigure = TRUE

# Read in data
r1  <- readRDS(file="results/overall_seroprevalence_and_R0_mab_rw_threshold_16_age_classes_40_m0_1_sample_size_500.rds")
r2  <- readRDS(file="results/overall_seroprevalence_and_R0_mab_rw_threshold_64_age_classes_40_m0_0.9_sample_size_500.rds")

r1$cutoff = "1:16"
r2$cutoff = "1:64"

r <- rbind(r1,r2)

q1 <- ggplot(data=r[r$Year %in% seq(2006,2017),], aes(x=Year)) +
  geom_line(aes(y=overall_seroprev_median, color=cutoff)) +
  geom_ribbon(aes(ymin=overall_seroprev_low,ymax=overall_seroprev_up,fill=cutoff)) +
  ylim(0,1) +
  scale_x_continuous(breaks=seq(2006,2017,by=2)) +  
  scale_color_manual(values=cutoff_cols) +
  scale_fill_manual(values=alpha(cutoff_cols,0.3)) +
  labs(y = "Overall seroprevalence") +
  theme_bw()  

q1 <- q1 + labs(color = "Cut-off",
                fill = "Cut-off") +
  theme(legend.position = c(.75, .35)) 

# Save
if(savefigure){
  ggsave(filename = "figures/Figure_4.pdf",q1,height=90,width=90,units="mm")
}



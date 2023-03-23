rm(list=ls())

library(cowplot)
library(ggpubr)

source("src/R/my_colors.R")

savefigure = TRUE


##-- Seroprevalence by age and year --##

# Read in data
p1  <- readRDS(file="figures/seroprevalence_by_age_and_year_mab_rw_threshold_16_age_classes_40_m0_1_sample_size_500.rds")
p2  <- readRDS(file="figures/seroprevalence_by_age_and_year_mab_rw_threshold_64_age_classes_40_m0_0.9_sample_size_500.rds")

# Add title
p1 <- p1 + labs(title="Seropositivity cut-off 1:16")
p2 <- p2 + labs(title="Seropositivity cut-off 1:64")


##-- Number of infections by age and year --##

# Read in data
q1  <- readRDS(file="figures/number_infections_by_age_and_year_mab_rw_threshold_16_age_classes_40_m0_1_sample_size_500.rds")
q2  <- readRDS(file="figures/number_infections_by_age_and_year_mab_rw_threshold_64_age_classes_40_m0_0.9_sample_size_500.rds")

q1 <- q1 + theme(legend.position="none")
q2 <- q2 + theme(legend.position="none")


##-- Mean age at infection --##

# Read in data
r1  <- readRDS(file="figures/mean_age_at_infection_mab_rw_threshold_16_age_classes_40_m0_1_sample_size_500.rds")
r2  <- readRDS(file="figures/mean_age_at_infection_mab_rw_threshold_64_age_classes_40_m0_0.9_sample_size_500.rds")

# Changes in x and y axis
r1 <- r1 + labs(y="Mean age",x="") + theme(axis.text.x = element_text(angle=90),
                                           plot.background = element_rect(fill='transparent', color=NA))
r2 <- r2 + labs(y="Mean age",x="") + theme(axis.text.x = element_text(angle=90),
                                           plot.background = element_rect(fill='transparent', color=NA))
# Insets
q1_with_r1 <-
  ggdraw() +
  draw_plot(q1) +
  draw_plot(r1, x = 0.55, y = 0.55, width = .4, height = .4)

q2_with_r2 <-
  ggdraw() +
  draw_plot(q2) +
  draw_plot(r2, x = 0.55, y = 0.55, width = .4, height = .4)


# Combine all
t <- ggarrange(p1, p2, q1_with_r1, q2_with_r2, ncol=2, nrow=2, common.legend = TRUE, legend="top",
               labels = c('A','B','C','D')) 
t

# Save
if(savefigure){
  ggsave(filename = "figures/Figure_5_panel.pdf",t,height=180,width=180,units="mm")
}


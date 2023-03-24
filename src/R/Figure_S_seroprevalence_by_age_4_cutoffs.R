rm(list=ls())

library(ggplot2)
library(cowplot)

savefigure = TRUE

# Read in data
# p_16 <- readRDS(file="figures/seroprevalence_all_ages_16.rds")
p_32 <- readRDS(file="figures/seroprevalence_all_ages_32.rds")
# p_64  <- readRDS(file="figures/seroprevalence_all_ages_64.rds")
p_128  <- readRDS(file="figures/seroprevalence_all_ages_128.rds")

# Add titles
# p1 <- p_16 + ggtitle("Seropositivity cut-off >=1:16") +
#   theme(legend.position = c(0.7, 0.4))
p2 <- p_32 + ggtitle("Seropositivity cut-off >=1:32") +
  # theme(legend.position = "none") 
  theme(legend.position = c(0.7, 0.4))
# p3 <- p_64 + ggtitle("Seropositivity cut-off >=1:64") +
#   theme(legend.position = "none")
p4 <- p_128 + ggtitle("Seropositivity cut-off >=1:128") +
  theme(legend.position = "none")

# Combine
p <- plot_grid(p2, p4, labels = c('A', 'B'), ncol=2, nrow=1)

# Save
if(savefigure){
  ggsave(filename = "figures/Figure_S_seroprevalence_2_cutoffs.pdf",p,height=90,width=180,units="mm")
}


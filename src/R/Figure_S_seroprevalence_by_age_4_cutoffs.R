rm(list=ls())

library(cowplot)

savefigure = TRUE

# Read in data
p_16 <- readRDS(file="figures/seroprevalence_all_ages_16.rds")
p_32 <- readRDS(file="figures/seroprevalence_all_ages_32.rds")
p_64  <- readRDS(file="figures/seroprevalence_all_ages_64.rds")
p_128  <- readRDS(file="figures/seroprevalence_all_ages_128.rds")

# Add titles
p1 <- p_16 + ggtitle("Seropositivity cut-off >=1:16") +
  theme(legend.position = c(0.7, 0.4))
p2 <- p_32 + ggtitle("Seropositivity cut-off >=1:32") +
  theme(legend.position = "none") 
p3 <- p_64 + ggtitle("Seropositivity cut-off >=1:64") +
  theme(legend.position = "none")
p4 <- p_128 + ggtitle("Seropositivity cut-off >=1:128") +
  theme(legend.position = "none")

# Combine
p <- plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), ncol=2, nrow=2)

# Save
if(savefigure){
  ggsave(filename = "figures/Figure_S_seroprevalence_4_cutoffs.pdf",p,height=180,width=180,units="mm")
}


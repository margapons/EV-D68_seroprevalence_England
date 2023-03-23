rm(list=ls())

library(cowplot)

source("src/R/my_colors.R")

savefigure = TRUE

# Read in data
p1  <- readRDS(file="figures/foi_mab_rw.rds")
p2  <- readRDS(file="figures/estimated_prop_w_mab_mab_rw.rds")

# Combine
p <- plot_grid(p1, p2, labels = c('A', 'B'), ncol=2, nrow=1)

p

# Save
if(savefigure){
  ggsave(filename = "figures/Figure_3_panel.pdf",p,height=90,width=180,units="mm")
}

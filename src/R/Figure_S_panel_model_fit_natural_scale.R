rm(list=ls())

library(cowplot)

savefigure = TRUE

model_abb = "mab_rw"

# Read in data
p1  <- readRDS(file=paste0("figures/model_fit_",model_abb,"_threshold_16_age_classes_40_m0_1.rds"))
p2  <- readRDS(file=paste0("figures/model_fit_",model_abb,"_threshold_64_age_classes_40_m0_0.9.rds"))

# Add titles
p1 <- p1 + ggtitle("Seropositivity cut-off 1:16") +
  scale_x_continuous(breaks=c(0,10,20,30,40),
                     labels=c(0,10,20,30,40)) 

p2 <- p2 + ggtitle("Seropositivity cut-off 1:64") +
  scale_x_continuous(breaks=c(0,10,20,30,40),
                     labels=c(0,10,20,30,40)) 

# Combine
p <- plot_grid(p1, p2, labels = c('A', 'B'), ncol=1, nrow=2)

p

# Save
if(savefigure){
  ggsave(filename = paste0("figures/Figure_model_fit_",model_abb,"_natural_scale.pdf"),p,height=130,width=183,units="mm")
}


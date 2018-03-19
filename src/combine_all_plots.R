library(cowplot)
#library(ggplot2)
library(dplyr)

source('plot_divergence_vs_time.R')
source('plot_divergence_vs_time_by_domain.R')

# Adjust individual plots
divergence_by_domain_plot <- divergence_by_domain_plot + 
  theme(strip.background = element_blank(), strip.text = element_blank())

glyc_vs_time_plot <- glyc_vs_time_plot + xlab('')

final_plot <- plot_grid(divergence_vs_time_plot, glyc_vs_time_plot, divergence_by_domain_plot,
          nrow = 3, labels = c('a','b','c'), label_y = c(1,1.1,1.1))



pdf('../figures/final_plot.pdf', height = 8.5, width = 10)
plot(final_plot)
dev.off()

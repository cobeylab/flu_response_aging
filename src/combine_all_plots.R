library(cowplot)
#library(ggplot2)
library(dplyr)

source('plot_divergence_vs_time.R')
source('plot_divergence_vs_time_by_domain.R')

# Adjust individual plots
divergence_by_domain_plot <- divergence_by_domain_plot + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12), strip.text = element_text(size = 12))

divergence_vs_time_plot <- divergence_vs_time_plot + theme(axis.title.y = element_text(size = 12), 
                                                           axis.text.y = element_text(size = 12),
                                                           strip.text = element_text(size = 12))

glyc_vs_time_plot <- glyc_vs_time_plot + xlab('') + theme(axis.title.y = element_text(size = 12), 
                                                           axis.text.y = element_text(size = 12),
                                                          strip.text = element_text(size = 12))

final_plot <- plot_grid(divergence_vs_time_plot, glyc_vs_time_plot, divergence_by_domain_plot,
          nrow = 3, labels = c('A','B','C'), label_y = c(1,1.1,1.1), label_size = 12)



pdf('../figures/final_plot.pdf', height = 8.5, width = 10)
plot(final_plot)
dev.off()

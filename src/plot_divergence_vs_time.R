library(cowplot)
#library(ggplot2)
library(dplyr)

# Import H3N2 results
H3N2_HA_genbank_results <- as_tibble(read.csv('../results/H3N2_HA_genbank_distances.csv', header= T))
H3N2_NA_gisaid_results <- as_tibble(read.csv('../results/H3N2_NA_gisaid_distances.csv', header= T))
H3N2_NP_genbank_results <- as_tibble(read.csv('../results/H3N2_NP_genbank_distances.csv', header= T))

# Merge tibbles
divergence_results <- bind_rows(mutate(H3N2_HA_genbank_results, subtype = 'H3N2', segment = 'HA', database = 'genbank'),
mutate(H3N2_NA_gisaid_results, subtype = 'H3N2', segment = 'NA', database = 'gisaid'),
mutate(H3N2_NP_genbank_results, subtype = 'H3N2', segment = 'NP', database = 'genbank'))

H3N2_plot <- ggplot(divergence_results, aes(x = year, y = distance_from_reference, 
                               color = segment, label = isolate_id)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(se = F) +
  xlab('Year') +
  ylab('Amino acid divergence from\nA/Hong Kong/1-5/1968') +
  scale_color_brewer(type = 'qual') +
  theme(legend.position = 'top')
  
# ggsave crashing for some reason
pdf('../figures/H3N2_divergence_vs_time.pdf', width = 8, height = 6)
plot(H3N2_plot)
dev.off()

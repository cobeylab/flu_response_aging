library(cowplot)
#library(ggplot2)
library(dplyr)

# Import H3N2 results (based on KC alignments)
#H3N2_HA_genbank_results <- as_tibble(read.csv('../results/H3N2_HA_genbank_distances.csv', header= T))
#H3N2_NA_gisaid_results <- as_tibble(read.csv('../results/H3N2_NA_gisaid_distances.csv', header= T))
#H3N2_NP_genbank_results <- as_tibble(read.csv('../results/H3N2_NP_genbank_distances.csv', header= T))

# Import H3N2 results 
H3N2_HA_genbank_results <- as_tibble(read.csv('../results/genbank_H3N2_HA_distances.csv', header= T))
H3N2_NA_genbank_results <- as_tibble(read.csv('../results/genbank_H3N2_NA_distances.csv', header= T))
H3N2_NP_genbank_results <- as_tibble(read.csv('../results/genbank_H3N2_NP_distances.csv', header= T))

# Import H1N1 results
H1N1_HA_genbank_results <- as_tibble(read.csv('../results/genbank_H1N1_HA_distances.csv', header= T))
H1N1_NA_genbank_results <- as_tibble(read.csv('../results/genbank_H1N1_NA_distances.csv', header= T))
H1N1_NP_genbank_results <- as_tibble(read.csv('../results/genbank_H1N1_NP_distances.csv', header= T))

# Filter H1N1 to keep only sequences from 1934 on (year of the reference sequence)
H1N1_HA_genbank_results <- filter(H1N1_HA_genbank_results, year >= 1977)
H1N1_NA_genbank_results <- filter(H1N1_NA_genbank_results, year >= 1977)
H1N1_NP_genbank_results <- filter(H1N1_NP_genbank_results, year >= 1977)

# Merge tibbles
divergence_results <- bind_rows(mutate(H3N2_HA_genbank_results, subtype = 'H3N2', segment = 'HA', database = 'genbank'),
                                mutate(H3N2_NA_genbank_results, subtype = 'H3N2', segment = 'NA', database = 'genbank'),
                                mutate(H3N2_NP_genbank_results, subtype = 'H3N2', segment = 'NP', database = 'genbank'),
                                mutate(H1N1_HA_genbank_results, subtype = 'H1N1', segment = 'HA', database = 'genbank'),
                                mutate(H1N1_NA_genbank_results, subtype = 'H1N1', segment = 'NA', database = 'genbank'),
                                mutate(H1N1_NP_genbank_results, subtype = 'H1N1', segment = 'NP', database = 'genbank')
                                ) %>% 
  mutate(subtype = factor(subtype, levels = c('H3N2','H1N1'))
         )

divergence_vs_time_plot <- ggplot(divergence_results, aes(x = year, y = distance_from_reference, 
                               color = segment, label = isolate_id)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(se = F) +
  facet_grid(.~subtype, scales = 'free_x') +
  xlab('Year') +
  ylab('Amino acid divergence from reference strain') +
  scale_color_brewer(type = 'qual') +
  theme(legend.position = 'top') 
  

# ggsave crashing for some reason
pdf('../figures/divergence_vs_time_plot.pdf', width = 10, height = 6)
plot(divergence_vs_time_plot)
dev.off()

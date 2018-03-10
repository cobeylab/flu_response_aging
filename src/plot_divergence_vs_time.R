library(cowplot)
#library(ggplot2)
library(dplyr)

# Import H3N2 results 
H3N2_NA_genbank_results <- as_tibble(read.csv('../results/genbank_H3N2_NA_distances.csv', header= T))
H3N2_NP_genbank_results <- as_tibble(read.csv('../results/genbank_H3N2_NP_distances.csv', header= T))
H3N2_HA_genbank_results <- as_tibble(read.csv('../results/genbank_H3N2_HA_distances.csv', header= T))

# Import H1N1 results
H1N1_HA_genbank_results <- as_tibble(read.csv('../results/genbank_H1N1_HA_distances.csv', header= T))
H1N1_NA_genbank_results <- as_tibble(read.csv('../results/genbank_H1N1_NA_distances.csv', header= T))
H1N1_NP_genbank_results <- as_tibble(read.csv('../results/genbank_H1N1_NP_distances.csv', header= T))

# Filter H1N1 to keep only sequences from 1934 on (year of the reference sequence)
H1N1_HA_genbank_results <- filter(H1N1_HA_genbank_results, year >= 1934)
H1N1_NA_genbank_results <- filter(H1N1_NA_genbank_results, year >= 1934)
H1N1_NP_genbank_results <- filter(H1N1_NP_genbank_results, year >= 1934)

# Import B Victoria results
B_Victoria_HA_genbank_results <- as_tibble(read.csv('../results/genbank_B_HA_Victoria_distances.csv', header= T))
B_Victoria_NA_genbank_results <- as_tibble(read.csv('../results/genbank_B_NA_Victoria_distances.csv', header= T))
B_Victoria_NP_genbank_results <- as_tibble(read.csv('../results/genbank_B_NP_Victoria_distances.csv', header= T))

# Filter B Victoria results to keep only sequences after 1987 (year of the reference strain)
B_Victoria_HA_genbank_results <- filter(B_Victoria_HA_genbank_results, year >= 1987)
B_Victoria_NA_genbank_results <- filter(B_Victoria_NA_genbank_results, year >= 1987)
B_Victoria_NP_genbank_results <- filter(B_Victoria_NP_genbank_results, year >= 1987)

# Import B Yamagata results
B_Yamagata_HA_genbank_results <- as_tibble(read.csv('../results/genbank_B_HA_Yamagata_distances.csv', header= T))
B_Yamagata_NA_genbank_results <- as_tibble(read.csv('../results/genbank_B_NA_Yamagata_distances.csv', header= T))
B_Yamagata_NP_genbank_results <- as_tibble(read.csv('../results/genbank_B_NP_Yamagata_distances.csv', header= T))

# Filter B Yamagata results to keep only sequences after 1988 (year of the reference strain)
B_Yamagata_HA_genbank_results <- filter(B_Yamagata_HA_genbank_results, year >= 1988)
B_Yamagata_NA_genbank_results <- filter(B_Yamagata_NA_genbank_results, year >= 1988)
B_Yamagata_NP_genbank_results <- filter(B_Yamagata_NP_genbank_results, year >= 1988)

# Merge tibbles
divergence_results <- bind_rows(mutate(H3N2_HA_genbank_results, subtype = 'H3N2', segment = 'HA', database = 'genbank'),
                                mutate(H3N2_NA_genbank_results, subtype = 'H3N2', segment = 'NA', database = 'genbank'),
                                mutate(H3N2_NP_genbank_results, subtype = 'H3N2', segment = 'NP', database = 'genbank'),
                                mutate(H1N1_HA_genbank_results, subtype = 'H1N1', segment = 'HA', database = 'genbank'),
                                mutate(H1N1_NA_genbank_results, subtype = 'H1N1', segment = 'NA', database = 'genbank'),
                                mutate(H1N1_NP_genbank_results, subtype = 'H1N1', segment = 'NP', database = 'genbank'),
                                mutate(B_Victoria_HA_genbank_results, subtype = 'B/Victoria', segment = 'HA', database = 'genbank'),
                                mutate(B_Victoria_NA_genbank_results, subtype = 'B/Victoria', segment = 'NA', database = 'genbank'),
                                mutate(B_Victoria_NP_genbank_results, subtype = 'B/Victoria', segment = 'NP', database = 'genbank'),
                                mutate(B_Yamagata_HA_genbank_results, subtype = 'B/Yamagata', segment = 'HA', database = 'genbank'),
                                mutate(B_Yamagata_NA_genbank_results, subtype = 'B/Yamagata', segment = 'NA', database = 'genbank'),
                                mutate(B_Yamagata_NP_genbank_results, subtype = 'B/Yamagata', segment = 'NP', database = 'genbank')
                                ) %>% 
  mutate(subtype = factor(subtype, levels = c('H3N2','H1N1', 'B/Victoria','B/Yamagata'))
         )

divergence_vs_time_plot <- ggplot(divergence_results, aes(x = year, y = distance_from_reference, 
                               color = segment, label = isolate_id)) + 
  geom_point(alpha = 0.5, size = 1) + 
  geom_smooth(se = F) +
  facet_grid(.~subtype, scales = 'free_x') +
  xlab('Year') +
  ylab('Amino acid divergence\nfrom reference strain') +
  scale_color_brewer(type = 'qual') +
  theme(legend.position = c(0.003,0.8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10)) 
  

# ggsave crashing for some reason
pdf('../figures/divergence_vs_time_plot.pdf', height = 3, width = 10)
plot(divergence_vs_time_plot)
dev.off()

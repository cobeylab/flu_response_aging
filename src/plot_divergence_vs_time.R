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

### AMINO ACID DIVERGENCE VS. TIME PLOT ###

divergence_vs_time_plot <- ggplot(divergence_results, aes(x = year, y = distance_from_reference)) + 
  geom_point(alpha = 0.5, size = 1.5, shape = 21, 
             aes(fill = segment, colour = segment), stroke = 0.1) + 
  #geom_smooth(se = F, aes(colour = segment)) +
  facet_grid(.~subtype, scales = 'free_x') +
  xlab('Year') +
  ylab('Amino acid divergence\nfrom reference strain') +
  scale_fill_brewer(type = 'qual') +
  scale_color_manual(values = c('springgreen4','darkorchid3','chocolate1')) +
theme(legend.position = c(0.003,0.78),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 10)) 

# Extract values from geom_point (to be able to ommit spline over gaps in the data)
divergence_spline_values <- ggplot_build(ggplot(divergence_results, aes(x = year, y = distance_from_reference)) +
  geom_smooth(se = F, aes(colour = segment)) + facet_grid(.~subtype))[[1]][[1]]

divergence_spline_values <- as_tibble(divergence_spline_values) %>%
  rename(year = x, spline_fit = y, subtype = PANEL, segment = group) %>%
  mutate(subtype = ifelse(subtype == '1', 'H3N2', subtype)) %>%
  mutate(subtype = ifelse(subtype == '2', 'H1N1', subtype)) %>%
  mutate(subtype = ifelse(subtype == '3', 'B/Victoria', subtype)) %>%
  mutate(subtype = ifelse(subtype == '4', 'B/Yamagata', subtype)) %>%
  mutate(subtype = factor(subtype, levels = c('H3N2','H1N1', 'B/Victoria','B/Yamagata'))) %>%
  mutate(segment = ifelse(segment == '1', 'HA', segment)) %>%
  mutate(segment = ifelse(segment == '2', 'NA', segment)) %>%
  mutate(segment = ifelse(segment == '3', 'NP', segment)) %>%
  mutate(segment = factor(segment, levels = c('HA','NA','NP')) )%>%
  mutate(spline_fit = ifelse(subtype == 'H1N1' & year >1957 & year <1977,
                             NA, spline_fit))

divergence_vs_time_plot <- divergence_vs_time_plot +
  geom_line(data = divergence_spline_values,aes(y = spline_fit,x = year,
                                                color = factor(segment)))
  
# ggsave crashing for some reason
pdf('../figures/divergence_vs_time_plot.pdf', height = 3, width = 10)
plot(divergence_vs_time_plot)
dev.off()

### PLOT WITH FRACTION OF GLYCOSYLATION SITES OVER TIME RELATIVE TO REFERENCE STRAIN ##
glyc_vs_time_plot <- ggplot(divergence_results, aes(x = year, y = fraction_glyc_seq - fraction_glyc_ref , 
                                                    color = segment)) + 
  geom_point(alpha = 0.5, size = 1.5, shape = 21, 
             aes(fill = segment, colour = segment), stroke = 0.1) + 
  #geom_smooth(se = F, aes(colour = segment)) +
  facet_grid(.~subtype, scales = 'free_x') +
  xlab('Year') +
  ylab('Change in fraction of\npotentially glycosylated sites') +
  scale_fill_brewer(type = 'qual') +
  scale_color_manual(values = c('springgreen4','darkorchid3','chocolate1')) +
  theme(legend.position = c(0.003,0.78),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10)) 

# Extract values from geom_point (to be able to ommit spline over gaps in the data)
glycosylation_spline_values <- ggplot_build(ggplot(divergence_results, aes(x = year, y = fraction_glyc_seq - fraction_glyc_ref)) +
                                           geom_smooth(se = F, aes(colour = segment)) + facet_grid(.~subtype))[[1]][[1]]

glycosylation_spline_values <- as_tibble(glycosylation_spline_values) %>%
  rename(year = x, spline_fit = y, subtype = PANEL, segment = group) %>%
  mutate(subtype = ifelse(subtype == '1', 'H3N2', subtype)) %>%
  mutate(subtype = ifelse(subtype == '2', 'H1N1', subtype)) %>%
  mutate(subtype = ifelse(subtype == '3', 'B/Victoria', subtype)) %>%
  mutate(subtype = ifelse(subtype == '4', 'B/Yamagata', subtype)) %>%
  mutate(subtype = factor(subtype, levels = c('H3N2','H1N1', 'B/Victoria','B/Yamagata'))) %>%
  mutate(segment = ifelse(segment == '1', 'HA', segment)) %>%
  mutate(segment = ifelse(segment == '2', 'NA', segment)) %>%
  mutate(segment = ifelse(segment == '3', 'NP', segment)) %>%
  mutate(segment = factor(segment, levels = c('HA','NA','NP')) )%>%
  mutate(spline_fit = ifelse(subtype == 'H1N1' & year >1957 & year <1977,
                             NA, spline_fit))

glyc_vs_time_plot <- glyc_vs_time_plot +
  geom_line(data = glycosylation_spline_values,aes(y = spline_fit,x = year,
                                                color = factor(segment)))


# ggsave crashing for some reason
pdf('../figures/new_glycosylation_sites_vs_time_plot.pdf', height = 3, width = 10)
plot(glyc_vs_time_plot)
dev.off()

library(cowplot)
#library(ggplot2)
library(dplyr)

divergence_results <- tibble()

for(subtype in c('H3N2','H1N1','Victoria','Yamagata')){
  for(domain in c('HA1','HA2','NAH','NAS','NAA')){
    
    if(subtype %in% c('H3N2','H1N1')){
      results <- as_tibble(read.csv(paste('../results/genbank', subtype, domain, 'distances.csv',sep ='_')))
    }else{
      results <- as_tibble(read.csv(paste('../results/genbank_B', domain, subtype, 'distances.csv',sep ='_')))
    }
    
    results <- mutate(results, subtype = ifelse(subtype %in% c('H1N1','H3N2'),subtype, paste('B/',subtype, sep = '')),
                      domain = domain)
    if(subtype == 'H1N1'){
      results <- filter(results, year >= 1934)
    }
    if(subtype == 'Victoria'){
      results <- filter(results, year >= 1987)
    }
    if(subtype == 'Yamagata'){
      results <- filter(results, year >= 1988)
    }
    divergence_results <- bind_rows(divergence_results, results)
  }
}
divergence_results <- mutate(divergence_results,
                             subtype = factor(subtype,levels = c('H3N2','H1N1','B/Victoria','B/Yamagata')))
divergence_results <- mutate(divergence_results, 
                             domain = factor(domain, levels = c('HA1','HA2','NAH','NAS','NAA')))

divergence_by_domain_plot <- ggplot(divergence_results, aes(x = year,
                                                            y = 100*distance_from_reference,
                                                            color = domain)) +
  geom_point(alpha = 0.5, size = 1.2, shape = 21, 
             aes(fill = domain, colour = domain), stroke = 0.1) +
  facet_grid(.~subtype, scales = 'free_x') +
  xlab('Year') +
  ylab('% amino acid divergence\nfrom reference strain') +
  #scale_fill_brewer(type = 'div', palette = 3)
  scale_fill_manual(values = c('#7fbf7b','#d9f0d3','#CC66FF','#e7d4e8','#99CCFF'),
                    labels = c('HA head','HA stem', 'NA head', 'NA non-head', 'NA active sites')) +
  scale_color_manual(values = c('#339900','#00CC66','#660099','#CC36FF','#3333FF'),
                     labels = c('HA head','HA stem', 'NA head', 'NA non-head', 'NA active sites')) +
  theme(legend.position = c(0.003,0.78),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        legend.key.height=unit(0.7,"line"))

# Calculate means over time
mean_divergence <- divergence_results %>% group_by(year,domain,subtype) %>%
  summarise(mean_divergence = mean(distance_from_reference, na.rm = T))

# Add mean divergence over time to plot
divergence_by_domain_plot <- divergence_by_domain_plot + geom_line(data = filter(mean_divergence, subtype != 'H1N1', subtype != 'B/Victoria'),
                                    aes(x = year, y = 100*mean_divergence, 
                                        colour = factor(domain))
                                    ) +
  geom_line(data = filter(mean_divergence, subtype == 'H1N1', year <= 1957),
            aes(x = year, y = 100*mean_divergence, 
                colour = factor(domain)), alpha = 0.7) +
  geom_line(data = filter(mean_divergence, subtype == 'H1N1', year >= 1976),
            aes(x = year, y = 100*mean_divergence, 
                colour = factor(domain)), alpha = 0.7) +
  geom_line(data = filter(mean_divergence, subtype == 'B/Victoria',grepl('NA', domain), year <= 2001),
            aes(x = year, y = 100*mean_divergence, 
                colour = factor(domain)), alpha = 0.7) +
  geom_line(data = filter(mean_divergence, subtype == 'B/Victoria', grepl('NA', domain), year > 2001),
            aes(x = year, y = 100*mean_divergence, 
                colour = factor(domain)), alpha = 0.7) +
  geom_line(data = filter(mean_divergence, subtype == 'B/Victoria', grepl('NA', domain) == F),
            aes(x = year, y = 100*mean_divergence, 
                colour = factor(domain)), alpha = 0.7)

# ggsave crashing for some reason
pdf('../figures/divergence_by_domain_plot.pdf', height = 3, width = 10)
plot(divergence_by_domain_plot)
dev.off()

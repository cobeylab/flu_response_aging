library(tidyverse)
library(cowplot)
library(stringr)

# Processes Blast output to select, for each isolate, the best matching strain (Yamagata/Victoria) based on the maximim bit score value from Blast

# Read output of blastn against Vic/87 and Yam/88 strains:
blast_output <- read.csv('../results/B_seq_blast.csv', header =T)
blast_output <- as.tibble(blast_output)

# Process output to produce a tible with the accession # of each strain and the most likely lineage (Yam or Vic) based on bit score
blast_output <- blast_output %>% group_by(query_id) %>% 
  filter(bit_score == max(bit_score)) %>% 
  select(query_id, subject_id) %>% 
  mutate(lineage = ifelse(grepl('Victoria',subject_id), 'Victoria', 'Yamagata')) %>%
  separate(query_id, into = c('accession_number','isolate_name','collection_date','country'), sep = '\\|') %>%
  select(isolate_name, accession_number, lineage)

# Write csv
write.csv(blast_output, file = '../results/B_lineage_assignments.csv', row.names = F)

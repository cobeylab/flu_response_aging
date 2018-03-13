#!/bin/bash

# HA1-HA2 partitions based on A/Puerto Rico/8-1/1934 Accn. ACV49545 
# Added 1 aa to ends of HA1 and HA2 to account for insertion in our alignment.
# Genbank annotation: HA1 = 18-343, HA2 = 344:565
ha1_sites=18-344
ha2_sites=345-566

# NAH (NA head) vs. NAS (NA non-head, inc. stem and transmembrane domains) partitions
# Within the head, sites in the enzymatic site are denoted NAA
# Based on Air 2011 compared with our alignment, head starting at conserved cysteine at position 92.
nah_sites=92-470
nas_sites=1-91
naa_sites=118,119,151,152,156,179,223,225,228,275,277,278,293,295,368,402,425


# Calculating divergence over time for H1N1 HA
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 all ../../results/alignments/genbank_H1N1_HA_alignment.fasta ../../results/genbank_H1N1_HA_distances.csv

# For NA
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 all ../../results/alignments/genbank_H1N1_NA_alignment.fasta ../../results/genbank_H1N1_NA_distances.csv

# For NP
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 all ../../results/alignments/genbank_H1N1_NP_alignment.fasta ../../results/genbank_H1N1_NP_distances.csv

# For HA (HA1)
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 $ha1_sites ../../results/alignments/genbank_H1N1_HA_alignment.fasta ../../results/genbank_H1N1_HA1_distances.csv

# For HA (HA2)
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 $ha2_sites ../../results/alignments/genbank_H1N1_HA_alignment.fasta ../../results/genbank_H1N1_HA2_distances.csv

# For NA (head)
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 $nah_sites ../../results/alignments/genbank_H1N1_NA_alignment.fasta ../../results/genbank_H1N1_NAH_distances.csv

# For NA (stem and transmembrane domains)
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 $nas_sites ../../results/alignments/genbank_H1N1_NA_alignment.fasta ../../results/genbank_H1N1_NAS_distances.csv

# For NA (enzymatic sites)
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 $naa_sites ../../results/alignments/genbank_H1N1_NA_alignment.fasta ../../results/genbank_H1N1_NAA_distances.csv


#!/bin/bash

# Calculating divergence over time for H1N1 HA
#python ../divergence_vs_time.py ../../results/alignments/genbank_H1N1_HA_alignment.fasta A/Puerto_Rico/8-1/1934 ../../results/genbank_H1N1_HA_distances.csv
python ../divergence_vs_time.py ../../results/alignments/genbank_H1N1_HA_alignment.fasta A/USSR/92/1977 ../../results/genbank_H1N1_HA_distances.csv

# For NA
#python ../divergence_vs_time.py ../../results/alignments/genbank_H1N1_NA_alignment.fasta A/Puerto_Rico/8-1/1934 ../../results/genbank_H1N1_NA_distances.csv
python ../divergence_vs_time.py ../../results/alignments/genbank_H1N1_NA_alignment.fasta A/USSR/92/1977 ../../results/genbank_H1N1_NA_distances.csv

# For NP
#python ../divergence_vs_time.py ../../results/alignments/genbank_H1N1_NP_alignment.fasta A/Puerto_Rico/8-1/1934 ../../results/genbank_H1N1_NP_distances.csv
python ../divergence_vs_time.py ../../results/alignments/genbank_H1N1_NP_alignment.fasta A/USSR/92/1977 ../../results/genbank_H1N1_NP_distances.csv



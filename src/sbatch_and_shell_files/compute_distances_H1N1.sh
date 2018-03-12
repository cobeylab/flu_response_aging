#!/bin/bash

# Calculating divergence over time for H1N1 HA
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 all ../../results/alignments/genbank_H1N1_HA_alignment.fasta ../../results/genbank_H1N1_HA_distances.csv

# For NA
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 all ../../results/alignments/genbank_H1N1_NA_alignment.fasta ../../results/genbank_H1N1_NA_distances.csv

# For NP
python ../divergence_vs_time.py A/Puerto_Rico/8-1/1934 all ../../results/alignments/genbank_H1N1_NP_alignment.fasta ../../results/genbank_H1N1_NP_distances.csv




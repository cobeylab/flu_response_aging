#!/bin/bash

# Calculating divergence over time for B HA Victoria
python ../divergence_vs_time.py B/Victoria/02/1987 all ../../results/alignments/genbank_B_HA_Victoria_alignment.fasta ../../results/genbank_B_HA_Victoria_distances.csv

# For NA
python ../divergence_vs_time.py B/Victoria/02/1987 all ../../results/alignments/genbank_B_NA_Victoria_alignment.fasta ../../results/genbank_B_NA_Victoria_distances.csv

# For NP
python ../divergence_vs_time.py B/Victoria/02/1987 all ../../results/alignments/genbank_B_NP_Victoria_alignment.fasta ../../results/genbank_B_NP_Victoria_distances.csv



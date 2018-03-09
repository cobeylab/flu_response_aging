#!/bin/bash

# Calculating divergence over time for H3N2 HA
python ../divergence_vs_time.py ../../results/alignments/genbank_H3N2_HA_alignment.fasta A/Hong_Kong/01/1968 ../../results/genbank_H3N2_HA_distances.csv

# For NA
python ../divergence_vs_time.py ../../results/alignments/genbank_H3N2_NA_alignment.fasta A/Hong_Kong/01/1968 ../../results/genbank_H3N2_NA_distances.csv

# For NP
python ../divergence_vs_time.py ../../results/alignments/genbank_H3N2_NP_alignment.fasta A/Hong_Kong/01/1968 ../../results/genbank_H3N2_NP_distances.csv



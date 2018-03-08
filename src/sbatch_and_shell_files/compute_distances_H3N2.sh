#!/bin/bash

# Calculating divergence over time for H3N2 HA
python ../divergence_vs_time.py ../../results/alignments/KC_alignments/H3N2_HA_genbank_alignment.fasta A/Hong_Kong/1-5/1968 ../../results/H3N2_HA_genbank_distances.csv

# For NA
python ../divergence_vs_time.py ../../results/alignments/KC_alignments/H3N2_NA_gisaid_alignment.fasta A/Hong_Kong/1-5/1968 ../../results/H3N2_NA_gisaid_distances.csv

# For NP
python ../divergence_vs_time.py ../../results/alignments/KC_alignments/H3N2_NP_genbank_alignment.fasta A/Hong_Kong/1-5/1968 ../../results/H3N2_NP_genbank_distances.csv


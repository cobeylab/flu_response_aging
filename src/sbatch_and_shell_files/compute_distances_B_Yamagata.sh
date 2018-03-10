#!/bin/bash

# Calculating divergence over time for B HA Yamagata
python ../divergence_vs_time.py ../../results/alignments/genbank_B_HA_Yamagata_alignment.fasta B/Yamagata/16/1988 ../../results/genbank_B_HA_Yamagata_distances.csv

# For NA
python ../divergence_vs_time.py ../../results/alignments/genbank_B_NA_Yamagata_alignment.fasta B/Yamagata/16/1988 ../../results/genbank_B_NA_Yamagata_distances.csv

# For NP
python ../divergence_vs_time.py ../../results/alignments/genbank_B_NP_Yamagata_alignment.fasta B/Yamagata/16/1988 ../../results/genbank_B_NP_Yamagata_distances.csv



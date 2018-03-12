#!/bin/bash

# Calculating divergence over time for B HA Yamagata
python ../divergence_vs_time.py B/Yamagata/16/1988 all ../../results/alignments/genbank_B_HA_Yamagata_alignment.fasta ../../results/genbank_B_HA_Yamagata_distances.csv

# For NA
python ../divergence_vs_time.py B/Yamagata/16/1988 all ../../results/alignments/genbank_B_NA_Yamagata_alignment.fasta ../../results/genbank_B_NA_Yamagata_distances.csv

# For NP
python ../divergence_vs_time.py B/Yamagata/16/1988 all ../../results/alignments/genbank_B_NP_Yamagata_alignment.fasta ../../results/genbank_B_NP_Yamagata_distances.csv



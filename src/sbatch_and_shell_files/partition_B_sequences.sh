#!/bin/bash

# Partitioning B HA sequences into Victoria and Yamagata, based on blast analysis of HA sequences
python ../partition_B_sequences.py ../../data/sequences/genbank_B_HA.fasta ../../results/B_lineage_assignments.csv

# Partitioning B NA sequences...
python ../partition_B_sequences.py ../../data/sequences/genbank_B_NA.fasta ../../results/B_lineage_assignments.csv

# Partitioning B Np sequences...
python ../partition_B_sequences.py ../../data/sequences/genbank_B_NP.fasta ../../results/B_lineage_assignments.csv


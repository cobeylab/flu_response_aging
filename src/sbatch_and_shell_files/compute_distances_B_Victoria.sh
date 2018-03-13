#!/bin/bash

# HA1 and HA2 sites from A/Brisbane/60/2008 PDB:4FQM
# Provided end for ha2 was 538, we set it to the last position in the alignment (585)
ha1_sites=16-362
ha2_sites=363-585

# NAH (NA head) vs. NAS (NA non-head, inc. stem and transmembrane domains) partitions
# Within the head, sites in the enzymatic site are denoted NAA
# Based on Air 2011 compared with our alignment visually for homology, head starting at conserved cysteine at position 88 instead of 92.
nah_sites=88-467
nas_sites=1-87
naa_sites=117,118,150,151,155,178,222,224,227,274,276,277,293,295,368,410,429


# Calculating divergence over time for B HA Victoria
python ../divergence_vs_time.py B/Victoria/02/1987 all ../../results/alignments/genbank_B_HA_Victoria_alignment.fasta ../../results/genbank_B_HA_Victoria_distances.csv

# For NA
python ../divergence_vs_time.py B/Victoria/02/1987 all ../../results/alignments/genbank_B_NA_Victoria_alignment.fasta ../../results/genbank_B_NA_Victoria_distances.csv

# For NP
python ../divergence_vs_time.py B/Victoria/02/1987 all ../../results/alignments/genbank_B_NP_Victoria_alignment.fasta ../../results/genbank_B_NP_Victoria_distances.csv

# HA (HA1)
python ../divergence_vs_time.py B/Victoria/02/1987 $ha1_sites ../../results/alignments/genbank_B_HA_Victoria_alignment.fasta ../../results/genbank_B_HA1_Victoria_distances.csv

# HA (HA2)
python ../divergence_vs_time.py B/Victoria/02/1987 $ha2_sites ../../results/alignments/genbank_B_HA_Victoria_alignment.fasta ../../results/genbank_B_HA2_Victoria_distances.csv

# NA (head)
python ../divergence_vs_time.py B/Victoria/02/1987 $nah_sites ../../results/alignments/genbank_B_NA_Victoria_alignment.fasta ../../results/genbank_B_NAH_Victoria_distances.csv

# NA (stem and transmembrane)
python ../divergence_vs_time.py B/Victoria/02/1987 $nas_sites ../../results/alignments/genbank_B_NA_Victoria_alignment.fasta ../../results/genbank_B_NAS_Victoria_distances.csv

# NA (enzymatic sites)
python ../divergence_vs_time.py B/Victoria/02/1987 $naa_sites ../../results/alignments/genbank_B_NA_Victoria_alignment.fasta ../../results/genbank_B_NAA_Victoria_distances.csv



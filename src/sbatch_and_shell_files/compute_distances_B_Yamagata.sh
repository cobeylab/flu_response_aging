#!/bin/bash

# HA1 and HA2 sites from A/Brisbane/60/2008 PDB:4FQM, a Yamagata strain
# Visually assessed homology to confirm that the sames positions apply to to Yam 
# Provided end for ha2 was 538, we set it to the last position in the alignment (585)
ha1_sites=16-362
ha2_sites=363-585

# NAH (NA head) vs. NAS (NA non-head, inc. stem and transmembrane domains) partitions
# Within the head, sites in the enzymatic site are denoted NAA
# Based on Air 2011 compared with our alignment visually for homology, head starting at conserved cysteine at position 87 instead of 92.
nah_sites=87-466
nas_sites=1-86
naa_sites=116,117,149,150,154,177,221,223,226,273,275,276,292,294,367,409,428

# Calculating divergence over time for B HA Yamagata
python ../divergence_vs_time.py B/Yamagata/16/1988 all ../../results/alignments/genbank_B_HA_Yamagata_alignment.fasta ../../results/genbank_B_HA_Yamagata_distances.csv

# For NA
python ../divergence_vs_time.py B/Yamagata/16/1988 all ../../results/alignments/genbank_B_NA_Yamagata_alignment.fasta ../../results/genbank_B_NA_Yamagata_distances.csv

# For NP
python ../divergence_vs_time.py B/Yamagata/16/1988 all ../../results/alignments/genbank_B_NP_Yamagata_alignment.fasta ../../results/genbank_B_NP_Yamagata_distances.csv

# HA (HA1)
python ../divergence_vs_time.py B/Yamagata/16/1988 $ha1_sites ../../results/alignments/genbank_B_HA_Yamagata_alignment.fasta ../../results/genbank_B_HA1_Yamagata_distances.csv

# HA (HA2)
python ../divergence_vs_time.py B/Yamagata/16/1988 $ha2_sites ../../results/alignments/genbank_B_HA_Yamagata_alignment.fasta ../../results/genbank_B_HA2_Yamagata_distances.csv

# NA (head)
python ../divergence_vs_time.py B/Yamagata/16/1988 $nah_sites ../../results/alignments/genbank_B_NA_Yamagata_alignment.fasta ../../results/genbank_B_NAH_Yamagata_distances.csv

# NA (stem and transmembrane)
python ../divergence_vs_time.py B/Yamagata/16/1988 $nas_sites ../../results/alignments/genbank_B_NA_Yamagata_alignment.fasta ../../results/genbank_B_NAS_Yamagata_distances.csv

# NA (enzymatic sites)
python ../divergence_vs_time.py B/Yamagata/16/1988 $naa_sites ../../results/alignments/genbank_B_NA_Yamagata_alignment.fasta ../../results/genbank_B_NAA_Yamagata_distances.csv

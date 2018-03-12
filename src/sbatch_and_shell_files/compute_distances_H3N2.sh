#!/bin/bash

# HA1-HA2 partitions based on A/Hong Kong/01/1968 Accn. AP039310 
ha1_sites=17-345
ha2_sites=346-566

# NAH (NA head) vs. NAS (NA non-head, inc. stem and transmembrane domains) partitions
# Within the head, sites in the enzymatic site are denoted NAA
# Based on Air 2011
nah_sites=92-469
nas_sites=1-92
naa_sites=118,119,151,152,156,178,222,224,227,274,276,277,292,294,371,406,425


# Calculating divergence over time for H3N2 HA (all sites)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 all ../../results/alignments/genbank_H3N2_HA_alignment.fasta ../../results/genbank_H3N2_HA_distances.csv

# For NA (all sites)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 all ../../results/alignments/genbank_H3N2_NA_alignment.fasta  ../../results/genbank_H3N2_NA_distances.csv

# For NP (all sites)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 all ../../results/alignments/genbank_H3N2_NP_alignment.fasta ../../results/genbank_H3N2_NP_distances.csv

# HA (HA1)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 $ha1_sites ../../results/alignments/genbank_H3N2_HA_alignment.fasta ../../results/genbank_H3N2_HA1_distances.csv

# HA (HA2)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 $ha2_sites ../../results/alignments/genbank_H3N2_HA_alignment.fasta ../../results/genbank_H3N2_HA2_distances.csv

# NA (head)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 $nah_sites ../../results/alignments/genbank_H3N2_NA_alignment.fasta  ../../results/genbank_H3N2_NAH_distances.csv

# NA (stem and transmembrane domains)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 $nas_sites ../../results/alignments/genbank_H3N2_NA_alignment.fasta  ../../results/genbank_H3N2_NAS_distances.csv

# NA (enzymatic sites)
python ../divergence_vs_time.py A/Hong_Kong/01/1968 $naa_sites ../../results/alignments/genbank_H3N2_NA_alignment.fasta  ../../results/genbank_H3N2_NAA_distances.csv


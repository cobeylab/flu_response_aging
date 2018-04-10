# flu_response_aging
Collaboration with Patrick's Wilson's group, led by Carole Henry: Elderly individuals no longer adapt their B cell responses to influenza efficiently

## Lineage assignment of influenza B sequences
Sbatch script ```blast_B_to_reference_seqs.sbatch``` contains the commands used for blasting influenza B HA  sequences in ```genbank_B_HA.fasta``` against the reference strains in ```reference_B_strains.fasta``` (both files are in ```data/sequences/```). Script ```process_blast_output.R``` processes the resulting csv file into another csv file listing the assigned lineage for each isolate (identified by its name and accession number). ```partition_B_sequences.py``` takes a fasta file with sequences   and the lineage assignment csv file produced by ```process_blast_output.R``` and outputs separate fasta files with sequences from the Yamagata and Victoria lineages. Shell script ```partition_B_sequences.sh``` executes ```partition_B_sequences.py``` for  ```genbank_B_HA.fasta```, ```genbank_B_NA.fasta``` and ```genbank_B_NP.fasta```.

## Sequence alignment
Commands for aligning HA, NA and NP sequences for each subtype / lineage using [MAFFT](https://mafft.cbrc.jp/alignment/software/) are contained in the ```*_alignment.sbatch``` sbatch files.

## Amino acid divergence and number of glycosylation sites
```divergence_vs_time.py``` uses functions implemented in ```alignment_analyses_functions.py``` to calculate % amino acid divergence and number of potential glycosylation sites over time given an alignment file, the name of a reference isolate and the range of sites to be considered. Shell scripts ```compute_distances*.sh``` execute that script for each subtype / lineage of influenza A and B.

## Plots
```combine_all_plots.R``` calls ```plot_divergence_vs_time.R``` and ```plot_divergence_vs_time_by_domain.R``` to produce the figure in the paper.

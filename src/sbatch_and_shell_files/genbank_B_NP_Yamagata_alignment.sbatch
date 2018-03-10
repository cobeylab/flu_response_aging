#!/bin/bash
#SBATCH --job-name=genbank_B_NP_Yamagata_alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --partition=cobey
#SBATCH -o genbank_B_NP_Yamagata_alignment.out        
#SBATCH -e genbank_B_NP_Yamagata_alignment.err      
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=20000

module load mafft/7.310

mafft --thread -8 --op 12 ../../data/sequences/genbank_B_NP_Yamagata.fasta > ../../results/alignments/genbank_B_NP_Yamagata_alignment.fasta

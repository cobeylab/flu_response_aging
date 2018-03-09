#!/bin/bash
#SBATCH --job-name=genbank_H3N2_NP_alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --partition=cobey
#SBATCH -o genbank_H3N2_NP_alignment.out        
#SBATCH -e genbank_H3N2_NP_alignment.err      
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=20000

module load mafft/7.310

mafft --op 3 --thread -8 ../../data/sequences/genbank_H3N2_NP.fasta > ../../results/alignments/genbank_H3N2_NP_alignment.fasta

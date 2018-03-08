#!/bin/bash
#SBATCH --job-name=genbank_H1N1_HA_alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --partition=cobey
#SBATCH -o genbank_H1N1_HA_alignment.out        
#SBATCH -e genbank_H1N1_HA_alignment.err      
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=20000

module load mafft/7.310

mafft --thread -8 ../../data/sequences/genbank_H1N1_HA.fasta > ../../results/alignments/genbank_H1N1_HA_alignment.fasta

#!/bin/bash
#SBATCH --job-name=gisaid_B_Victoria_HA_alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --partition=cobey
#SBATCH -o gisaid_B_Victoria_HA_alignment.out        
#SBATCH -e gisaid_B_Victoria_HA_alignment.err      
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=20000

module load mafft/7.310

mafft --op 5 --thread -8 ../../data/sequences/gisaid_B_Victoria_HA.fasta > ../../results/alignments/gisaid_B_Victoria_HA_alignment.fasta

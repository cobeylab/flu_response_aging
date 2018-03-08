#!/bin/bash
#SBATCH --job-name=gisaid_H3N2_HA_alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --partition=cobey
#SBATCH -o gisaid_H3N2_HA_alignment.out        
#SBATCH -e gisaid_H3N2_HA_alignment.err      
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=20000

module load mafft/7.310

mafft --thread -8 --op 5 ../../data/sequences/gisaid_H3N2_HA.fasta > ../../results/alignments/gisaid_H3N2_HA_alignment.fasta

#!/usr/bin/python
"""Partitions input fasta file with influenza B sequences from any segment into separate fasta files for Victoria and Yamagata lineages. Does so based on a csv file where isolate names are listed with their corresponding lineage based on blast"""

import sys
from Bio import SeqIO
import csv


# Read blast lineage assignment as a dictionary
lineage_assignment_path = '../results/B_lineage_assignments.csv'
lineage = {}
with open(lineage_assignment_path, 'r') as lineage_assignment_file:
    lineage_assignment = csv.reader(lineage_assignment_file)

    for line in lineage_assignment:
        isolate_name = line[0]
        isolate_lineage = line[2]
        lineage[isolate_name] = isolate_lineage

def main(argv):

    fasta_file_path = str(argv[1])

    output_path_vic = fasta_file_path.replace('.fasta', '_Victoria.fasta')
    output_path_yam = fasta_file_path.replace('.fasta', '_Yamagata.fasta')

    # Read and parse original fasta file
    with open(fasta_file_path, 'rU') as fasta_file:
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))

    victoria_sequences = []
    yamagata_sequences = []

    # For each sequence in fasta file...
    for record in sequences:
        # Extract isolate name
        isolate_name = record.id.split('|')[1]

        if isolate_name in lineage.keys():
            # Find isolate lineage from lineage dictionary
            isolate_lineage = lineage[isolate_name]

            # Store record in corresponding list for Victoria or Yamagata
            if isolate_lineage == 'Victoria':
                victoria_sequences.append(record)
            elif isolate_lineage == 'Yamagata':
                yamagata_sequences.append(record)

        else:
            print 'Warning: isolate ' + isolate_name + ' not present in lineage assignment csv file.'

    # Output Yam and Vic lists of sequences to different files:
    SeqIO.write(victoria_sequences, output_path_vic, 'fasta')
    SeqIO.write(yamagata_sequences, output_path_yam, 'fasta')

if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)



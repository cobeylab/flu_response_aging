#!/usr/bin/python
"""Partitions input fasta file with influenza B sequences from any segment into separate fasta files for Victoria and Yamagata lineages. Does so based on a csv file where isolate names are listed with their corresponding lineage based on blast"""

import sys
from Bio import SeqIO
import csv


def main(argv):

    fasta_file_path = str(argv[1])
    lineage_assignment_path = str(argv[2])

    # Read blast lineage assignment as a dictionary
    lineage = {}
    with open(lineage_assignment_path, 'r') as lineage_assignment_file:
        lineage_assignment = csv.reader(lineage_assignment_file)

        for line in lineage_assignment:
            isolate_name = line[0]
            isolate_lineage = line[2]
            lineage[isolate_name] = isolate_lineage

    output_path_vic = fasta_file_path.replace('.fasta', '_Victoria.fasta')
    output_path_yam = fasta_file_path.replace('.fasta', '_Yamagata.fasta')
    output_path_unassigned = fasta_file_path.replace('.fasta', '_unassigned.txt')

    # Read and parse original fasta file
    with open(fasta_file_path, 'rU') as fasta_file:
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))

    n_sequences = len(sequences)

    victoria_sequences = []
    yamagata_sequences = []
    unassigned_sequences = []

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
            unassigned_sequences.append(isolate_name)

    # Output Yam and Vic lists of sequences to different files:
    SeqIO.write(victoria_sequences, output_path_vic, 'fasta')
    SeqIO.write(yamagata_sequences, output_path_yam, 'fasta')

    with open(output_path_unassigned, 'w') as output_file_unassigned:
        output_file_unassigned.write(str(len(unassigned_sequences)) + ' of ' + str(n_sequences))
        output_file_unassigned.write(' (' + str(100 * float(len(unassigned_sequences)) / n_sequences) + '%) ')
        output_file_unassigned.write('isolates not present in the lineage assignment file. They are listed below.\n')

        for name in unassigned_sequences:
            output_file_unassigned.write(name + '\n')


if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)



#!/usr/bin/python
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
import sys

from alignment_analyses_functions import extract_year, extract_reference_seq, distance_from_reference, has_frameshift
from alignment_analyses_functions import fraction_sites_glycosylated


def main(argv):

    # Name of the reference isolate
    ref_isolate_name = str(argv[1])

    # Sites to be included in divergence calculation. Either "all" or a range in the format "X-Y"
    sites = str(argv[2])

    # Path to alignment file
    alignment_file_path = str(argv[3])

    # Path to output file
    output_file_path = str(argv[4])

    # Read alignment
    with open(alignment_file_path, 'r') as alignment_file, open(output_file_path, 'w') as output_file:

        # Write header to CSV file
        output_file.write('isolate_id,year,distance_from_reference,n_glyc_ref,n_glyc_seq,')
        output_file.write('fraction_glyc_ref, fraction_glyc_seq\n')

        alignment = AlignIO.read(alignment_file, 'fasta')

        ref_seq = extract_reference_seq(ref_isolate_name, alignment)

        if sites == 'all':
            sites = range(len(ref_seq))
        elif '-' in sites:
            sites = sites.split('-')
            sites = [int(site) - 1 for site in sites]
            sites = range(sites[0], sites[1] + 1)
        elif ',' in sites:
            sites = sites.split(',')
            sites = [int(site) - 1 for site in sites]

        # Fraction of sites potentially glycosylated in reference sequence:
        n_glyc_ref = str(fraction_sites_glycosylated(ref_seq)['n_glyc'])
        fraction_glyc_ref = str(fraction_sites_glycosylated(ref_seq)['fraction_glyc'])

        # Translate nucleotide sequences into amino acid sequences
        for record in alignment:
            # Replace all 'Xs' in the sequence with 'N'
            record.seq = Seq(str(record.seq).replace('x','n'), SingleLetterAlphabet())
            seq = record.seq

            seq_id = record.id
            seq_year = extract_year(seq_id)

            if not has_frameshift(str(seq)):
                seq = str(seq.translate(gap = '-'))
                seq_distance_from_ref = str(distance_from_reference(seq = seq, ref_seq = ref_seq, sites = sites))

                # Fraction of sites glycosylated in sequence
                n_glyc_seq = str(fraction_sites_glycosylated(seq)['n_glyc'])
                fraction_glyc_seq = str(fraction_sites_glycosylated(seq)['fraction_glyc'])

            else:
                seq_distance_from_ref = 'NA'
                fraction_glyc_seq = 'NA'
                n_glyc_seq = 'NA'

            output_file.write(','.join([seq_id, seq_year, seq_distance_from_ref,n_glyc_ref, n_glyc_seq,
                                        fraction_glyc_ref, fraction_glyc_seq]))
            output_file.write('\n')

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)











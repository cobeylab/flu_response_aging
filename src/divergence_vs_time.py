#!/usr/bin/python
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
import re
import sys

def extract_year(seq_id):
    if seq_id.find('unknown') == -1:
        date = re.search(r'\|[0-9]*\/[0-9]*\/[0-9]*', seq_id).group()
        date = date.replace('|','')
        year = date.split('/')[0]
        if len(year) == 0:
            year = 'NA'
        else:
            assert len(year) == 4
    else:
        year = 'NA'
    return (year)

def has_frameshift(seq):
    frameshift = False
    for codon_position in range(len(seq) / 3):
        codon = seq[codon_position * 3: codon_position * 3 + 3]
        n_gaps = len([base for base in codon if base == '-'])

        # Sequence has a frameshift mutation if any codon has 1 or 2 gaps (but not 3)
        if n_gaps == 1 or n_gaps == 2:
            frameshift = True

    return frameshift

def extract_reference_seq(ref_isolate_name, alignment):
    # Find index of SeqRecord that contains the reference sequence id in its label
    matching_records = []
    for record in alignment:
        isolate_name = record.id.split('|')[1]
        if isolate_name == ref_isolate_name:
            matching_records.append(record)
    assert len(matching_records) != 0, 'Reference isolate not found in alignment'

    assert len(matching_records) < 3, 'More than two occurrences of reference isolate found in alignment'

    # If reference strain is duplicated in the alignment, check that the sequences are identical at the A.A. level
    if len(matching_records) == 2:
        ref_record_1 = matching_records[0]
        ref_record_2 = matching_records[1]

        # Convert gaps to 'n' to avoid conflict with biopython -- will be translated as X and ignored anyway
        ref_seq_1 = Seq(str(ref_record_1.seq).replace('-','n'), SingleLetterAlphabet())
        ref_seq_2 = Seq(str(ref_record_2.seq).replace('-','n'), SingleLetterAlphabet())

        assert has_frameshift(str(ref_seq_1)) == False and has_frameshift(str(ref_seq_2)) == False, 'Duplicate reference sequences have a frameshift mutation'

        ref_seq_1 = str(ref_seq_1.translate())
        ref_seq_2 = str(ref_seq_2.translate())

        n_diffs = len([i for i in range(len(ref_seq_1)) if ref_seq_1[i] != ref_seq_2[i]])

        ref_seq = ref_seq_1

        if n_diffs >0:
            print 'Warning: there are duplicate reference sequences differing by ' + str(n_diffs) + ' amino acids. Choosing first listed.'

    elif len(matching_records) == 1:
        ref_record = matching_records[0]
        # Convert gaps to 'n' to avoid conflict with biopython -- will be translated as X and ignored anyway
        ref_seq = Seq(str(ref_record.seq).replace('-','n'), SingleLetterAlphabet())
        assert has_frameshift(str(ref_seq)) == False, 'Reference sequence has a frameshift mutation'
        ref_seq = str(ref_seq.translate())

    return ref_seq


def distance_from_reference(seq, ref_seq):
    # Disregards sites with gaps
    assert len(seq) == len(ref_seq)

    # Number of sites with gaps or ambiguities in one sequence or the other
    invalid_sites = [i for i in range(len(seq)) if seq[i] in {'-','X'} or ref_seq[i] in {'-','X'}]

    diff_sites = [i for i in range(len(seq)) if seq[i] != ref_seq[i]]

    diff_sites = [site for site in diff_sites if site not in invalid_sites]


    distance = float(len(diff_sites)) / (len(seq) - len(invalid_sites))
    return distance

def main(argv):

    alignment_file_path = str(argv[1])
    ref_isolate_name = str(argv[2])
    output_file_path = str(argv[3])

    # Read alignment
    with open(alignment_file_path, 'r') as alignment_file, open(output_file_path, 'w') as output_file:

        # Write header to CSV file
        output_file.write('isolate_id,year,distance_from_reference\n')

        alignment = AlignIO.read(alignment_file, 'fasta')


        ref_seq = extract_reference_seq(ref_isolate_name, alignment)

        # Translate nucleotide sequences into amino acid sequences
        for record in alignment:
            # Replace all 'Xs' in the sequence with 'N'
            record.seq = Seq(str(record.seq).replace('x','n'), SingleLetterAlphabet())
            seq = record.seq

            seq_id = record.id
            seq_year = extract_year(seq_id)

            if has_frameshift(str(seq)) == False:
                seq = str(seq.translate(gap = '-'))
                seq_distance_from_ref = str(distance_from_reference(seq = seq, ref_seq = ref_seq))
            else:
                seq_distance_from_ref = 'NA'

            output_file.write(','.join([seq_id, seq_year, seq_distance_from_ref]))
            output_file.write('\n')

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)











from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
import re

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


def distance_from_reference(seq, ref_seq, sites):
    """
    :param seq: string with query amino acid sequence
    :param ref_seq: string with reference amino acid sequence
    :param sites: list of sites (in python indexing) to be considered
    :return: amino acid divergence between query and ref. seqs. (n diffs. / number of sites excluding gaps and ambiguities)
    """

    assert len(seq) == len(ref_seq)

    # Number of sites with gaps or ambiguities in one sequence or the other
    invalid_sites = [i for i in sites if seq[i] in {'-','X'} or ref_seq[i] in {'-','X'}]

    # Number of sites that are different between sequence and reference sequence
    diff_sites = [i for i in sites if seq[i] != ref_seq[i]]

    # Exclude invalid sites (with gaps or ambiguities)
    diff_sites = [site for site in diff_sites if site not in invalid_sites]

    # Normalize difference by length of sequence excluding invalid sites
    distance = float(len(diff_sites)) / (len(sites) - len(invalid_sites))

    return distance

assert(distance_from_reference('MAANRSCNASCNPSC','MXANPSCDAS-NPSA', sites = range(15)) == float(3)/13)
assert(distance_from_reference('MAANRSCNASCNPSC','MXANPSCDAS-NPSA', sites = [0,1,4]) == float(1)/2)

def fraction_sites_glycosylated(seq):
    """
    :param seq: an amino acid sequence
    :return: number of glycosylation motifs, N-X-[ST]-X, where X is not proline
    """
    glycount = 0

    # Sites not to be included in the denominator of n glyc. sites per site
    n_invalid_sites = 0
    for i in range(len(seq) - 3):

        motif = seq[i:i+4]

        if motif[0] == 'N':
            if motif[1] != 'P' and motif[2] in ['S','T'] and motif[3] != 'P':
                glycount += 1

        # If site is an ambiguous amino acid (X) or a gap, exclude from denominator of # glyc. sites per site
        elif motif[0] == 'X' or motif[0] == '-':
            n_invalid_sites += 1

    # Fraction of sites glycosylated (denominator excludes gaps and ambiguous amino acids)
    fraction_glyc = float(glycount) / (len(seq) - n_invalid_sites)

    return {'n_glyc': glycount, 'fraction_glyc': fraction_glyc}

assert(fraction_sites_glycosylated('MAANRSCNASCNPSC')['n_glyc'] == 2)
assert(fraction_sites_glycosylated('MAANRSCNASCNPSC')['fraction_glyc'] == float(2) / 15)

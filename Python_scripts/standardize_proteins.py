"""
A script to standardize proteins from sys.argv[1].
Converts non-standard characters to 'X', discards proteins with less than 30 or more than 2500 aa,
converts all characters to uppercase.
Saves the result to sys.argv[2].
"""

from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACProtein
from matplotlib import pyplot as plt
import numpy as np
import copy, sys

def standardize_protein_sequences(sequences, alphabet=IUPACProtein.letters):
    """
    Returns a list of sequences in which all the residues not found in the supplied alphabet are changed to 'X',
    regardless whether they're lower- or upper-case.
    Note that BLAST denotes unknown nucleotides as 'N' instead of 'X', so this function should only be used for proteins.
    Gaps are left unchanged even if they're not found in the alphabet.
    :param sequences: List of SeqRecord objects.
    :param alphabet: A string consisting of the accepted letters.
    :return: List of modified SeqRecord objects.
    """
    sequences = copy.deepcopy(sequences)
    alphabet = str(alphabet).upper()
    for seqrec_nb, seqrec in enumerate(sequences):
        for i, c in enumerate(seqrec.seq):
            c = c.upper()  
            if c != '-' and c not in alphabet:
                newseq = seqrec.seq.tomutable()
                newseq[i] = 'X'
                newseq = newseq.toseq()
                seqrec.seq = newseq
    return sequences

LENGTH_RANGE = (30, 2500)   # range of acceptable protein length

infile, outfile = sys.argv[1:3]

all_prots = list(SeqIO.parse(infile, 'fasta'))
print('Loaded', len(all_prots), 'proteins')
all_prots = standardize_protein_sequences(all_prots)
# print(all_prots)
lengths = np.array([len(x) for x in all_prots])

aa_nbs = np.array([sum(c != 'X' and c.isupper() for c in x) for x in all_prots])
# print(lengths, aa_nbs)
#plt.subplot(121)
#plt.hist(np.log10(lengths), bins=100)
#plt.subplot(122)
#plt.hist(np.log10(aa_nbs[aa_nbs>0]), bins=100)
#plt.show()

#plt.hist(aa_nbs/lengths, bins=100)
#plt.show()

too_short = aa_nbs < 30
too_long = lengths > 2000
low_complexity = aa_nbs/lengths <= 0.5
print(sum(too_short), 'too short')
print(sum(too_long), 'too long')
print(sum(low_complexity), 'low complexity')
##1102 too short
##5011 too long
##739 low complexity
discard = too_short + too_long + low_complexity
retained = [p for p, d in zip(all_prots, discard) if not d]
SeqIO.write(retained, outfile,  'fasta')
    
    

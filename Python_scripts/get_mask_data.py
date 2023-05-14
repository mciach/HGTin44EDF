from Bio import SeqIO
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

def get_mask_info(seq):
    aa_nb = len(seq)
    lc_nb = sum(c.islower() for c in seq)
    lc_prop = lc_nb / aa_nb
    return (aa_nb, lc_nb, lc_prop)

S = SeqIO.parse(infile, 'fasta')
lc_data = [['ID', 'Length', 'Low complexity number', 'Low complexity proportions']]
for s in S:
    lc = get_mask_info(s)
    lc_data.append((s.id, ) + tuple(map(str, lc)))

with open(outfile, 'w') as h:
    for l in lc_data:
        h.write('\t'.join(l)+'\n')


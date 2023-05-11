import os, re

gff_files = ['GFFs/'+f for f in os.listdir('GFFs') if 'gff' == f[-3:]]

acc2intron_nb = {}
for f in gff_files:
    print('Processing', f)
    with open(f) as h:
        for l in h:
            if l[0] == '#':
                continue
            l = l.strip().split('\t')
            if l[2] == 'CDS':
                acc = re.search('protein_id=([^;]*)(;|$)', l[-1])
                if acc:
                    acc = acc.group(1)
                    assert ';' not in acc
                    assert '=' not in acc
                    assert '.' in acc
                    try:
                        acc2intron_nb[acc] += 1
                    except KeyError:
                        acc2intron_nb[acc] = 0
print('Verifying consistency')
basal_acc = [l.strip() for l in open('basal_accessions.txt') if l.strip()]
assert all(acc in acc2intron_nb for acc in basal_acc)

with open('intron_count.tsv' ,'w') as h:
    for key in acc2intron_nb:
        h.write(key + '\t' + str(acc2intron_nb[key]) + '\n')
        

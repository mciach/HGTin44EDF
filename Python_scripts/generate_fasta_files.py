"""
Skrypt do tworzenia plików FASTA z sekwencjami xenologów wszystkich oraz grzybowych.
Bierze tylko displaced_clade z warunkiem repositionable == Yes.
"""

from Bio import SeqIO
from ete3 import NCBITaxa, Tree
from path import Path
ncbi = NCBITaxa()

# inputy:
results_dir = Path('Results')

transferred_clusters = results_dir / 'hgt_results.tsv'
tree_directory = 'third_round_trees/'
blasted_proteins = 'standardized_final_sequences.fa'
blasted_taxid_table = 'final_homolog_taxid_table'

# outputy:
result_all = results_dir / 'all_hgt_sequences.fa'
##result_fungi = results_dir / 'fungal_hgt_sequences.fa'

# processing
fungi_taxid = '4751'

# bierzemy accessions wszystkich bialek  z rodzin xenologicznych
# uwzgledniamy bialka bakteryjne jako reference groups
accessions = []
with open(transferred_clusters) as h:
    for l in h:
        l = l.strip().split('\t')
        if not l or l[0] == 'cluster_id': continue
        cluster_id, criterion, well_supported = l[:3]
        if well_supported != 'Yes':
            continue
        T = Tree(tree_directory + cluster_id + '.treefile', format=1)
        for leaf in T:
            # assert leaf.name not in accessions
            accessions.append(leaf.name)
accessions = set(accessions)
# bierzemy sekwencje bialek z rodzin xenologicznych
proteins = []
S = SeqIO.parse(blasted_proteins, 'fasta')
for s in S:
    if s.id in accessions:
        proteins.append(s)
assert len(proteins) == len(accessions)

# zapisujemy wszystkie rodziny xenologiczne
SeqIO.write(proteins, result_all, 'fasta')

### Wybieramy sekwencje xenologow grzybowych
##fungal_accessions = set()
##with open(blasted_taxid_table) as h:
##    for l in h:
##        l = l.strip().split('\t')
##        if not l: continue
##        if l[0] in accessions:
##            species_taxid = l[1]
##            lineage = [str(lin) for lin in ncbi.get_lineage(species_taxid)]
##            if fungi_taxid in lineage:
##                fungal_accessions.add(l[0])
##fungal_sequences = [s for s in proteins if s.id in fungal_accessions]
##SeqIO.write(fungal_sequences, result_fungi, 'fasta')
##

# Laczymy sekwencje i accessions
##all_proteins = fungal_proteins
##all_proteins.extend(additional_proteins)  # uwaga - to modyfikuje fungal_proteins
##                                                                           # (chcemy oszczedzic pamiec)
##all_accessions = fungal_protein_accessions | additional_accessions
##assert len(all_accessions) == len(all_proteins)

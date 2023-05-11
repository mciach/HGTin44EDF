"""
A script to remove homologs with difficult taxonomy from global blast results.
Removes all homologs that do not belong to eukarya, archaea or bacteria,
also all non-fungal homologs with "incertae sedis" or "unclassified" anywhere in their lineage.
Input = table with accession and taxid.
Output = table with accession, corresponding species taxid, and the original taxid. 
Reads input from sys.argv[1] and saves the result in sys.argv[2].
"""
import sys
from ete3 import NCBITaxa
ncbi = NCBITaxa()

cellular_org_taxid = '131567'
fungi_taxid = '4751'

with open(sys.argv[1]) as h:
    original_taxid_set = set(l.strip().split()[1] for l in h)
print(len(original_taxid_set), 'unique taxIDs read')

try:
    outfile = sys.argv[2]
except IndexError:
    print('Provide the output filename!')
    raise

### Tree-based version:

##taxid_tree = ncbi.get_topology(taxid_set)
##cellular_orgs = taxid_tree.search_nodes(name=cellular_org_taxid)
##assert len(cellular_orgs) == 1
##cellular_orgs = cellular_orgs[0]
##print(len(cellular_orgs), 'taxIDs corresponding to cellular organisms')
##
##fungi = cellular_orgs.search_nodes(name=fungi_taxid)
##assert len(fungi) == 1
##fungi = fungi[0]
##print(len(fungi), 'taxIDs corresponding to fungi')
##
##leaf_fn = lambda x: x.is_leaf() or x == fungi or 'incertae' in x.sci_name.lower() or 'unclassified' in x.sci_name.lower()
##for n in cellular_orgs.traverse(strategy='postorder', is_leaf_fn = leaf_fn):
##    if 'incertae' in n.sci_name.lower() or 'unclassified' in n.sci_name.lower():
##        n.detach()
##        print('detached', n.sci_name, 'with', len(n), 'leaf taxIDs')
##        
##print(len(cellular_orgs), 'taxIDs remaining')
##fungi = cellular_orgs.search_nodes(name=fungi_taxid)
##assert len(fungi) == 1
##fungi = fungi[0]
##print(len(fungi), 'taxIDs corresponding to fungi remaining')

### Lineage-based version:
lineages = ncbi.get_lineage_translator(original_taxid_set)
print('Obtained', len(lineages), 'lineages')
# remove non-cellular organisms:
lineages = {tx: lineages[tx] for tx in lineages if lineages[tx][1] == 131567}
print('Retained', len(lineages), ' lineages of cellular organisms')
# remove 'unclassified' and 'incertae sedis' of non-fungal organisms:
taxid_list = [tx for tx in lineages]
species_taxids = {}
to_remove = [False]*len(taxid_list)
no_species = 0
fungi = 0
incertae = 0
unclassified = 0
environmental = 0
uncultured = 0
for i, tx in enumerate(taxid_list):
    is_fungus = False
    lng = lineages[tx]
    ranks = ncbi.get_rank(lng)
    ranks = [ranks[x] for x in lng]
    if not 'species' in ranks:
        to_remove[i] = True
        no_species += 1
        continue
    if lng[4] == 4751:  # fungus
        # we don't reject fungi incertae sedis
        fungi += 1
        is_fungus = True
    species_taxids[str(tx)] = lng[ranks.index('species')]
    names = ncbi.translate_to_names(lng)
    for n in names:
        if 'uncultured' in n:
            print(n)
        n = n.lower()
        if 'unclassified' in n:
            unclassified += 1
            to_remove[i] = True
            break
        elif 'environmental' in n:
            environmental += 1
            to_remove[i] = True
            break
        elif 'uncultured' in n:
            uncultured += 1
            to_remove[i] = True
            break
        elif not is_fungus and 'incertae' in n:
            incertae += 1
            to_remove[i] = True
            break
   
taxid_list = {str(tx) for tx, rm in zip(taxid_list, to_remove) if not rm}
print('Removed', no_species, 'lineages with no species')
print('Detected', fungi, 'fungal proteins')
print('Removed', incertae, 'incertae sedis lineages')
print('Removed', unclassified, 'unclassified lineages')
print('Removed', environmental, 'environmental lineages')
print('Removed', uncultured, 'uncultured lineages')
print('Retained', len(taxid_list), 'taxIDs with proper taxonomy,')
print('corresponding to', len({species_taxids[tx] for tx in taxid_list}), 'species')
# tu by jeszcze trzeba odfiltrować dziwne grzyby...

remaining_accessions = set()
nb_of_all_accessions = 0
with open(sys.argv[1]) as h:
    for l in h:
        nb_of_all_accessions += 1
        acc, txid = l.strip().split()
        if txid in taxid_list:
            remaining_accessions.add(str(acc) + '\t' + str(species_taxids[txid]) + '\t' + str(txid))
with open(outfile, 'w') as h:
    h.write('\n'.join(remaining_accessions) + '\n')
print(len(remaining_accessions), 'out of', nb_of_all_accessions, 'retained and saved to', outfile)

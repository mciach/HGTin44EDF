from Bio import SeqIO
import os, path

accessions = open('single_family_accessions.txt')
accessions = {l.strip() for l in accessions}

proteins = SeqIO.parse('standardized_proteins.fa', 'fasta')
sequences = []
for p in proteins:
    if p.id in accessions:
        accessions.remove(p.id)
        sequences.append(p)

SeqIO.write(sequences, 'global_blast_query.fa', 'fasta')




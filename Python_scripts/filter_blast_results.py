"""
Process blast results in clustering_blast_results.
Take only HSPs with at least 50% query cover and save in query_filtered_blast_results.
Take only HSPs with at least 50% query and subject cover and save in twoway_filtered_blast_results.
"""

query_length = open('sequence_length_table.txt')
query_length = {l.strip().split()[0]: int(l.strip().split()[1]) for l in query_length if l}
low_identity_discarded = 0
query_filter_discarded = 0
twoway_filter_discarded = 0
print('Sequence lengths read') 
blast_results = open('clustering_blast_results')
query_filter_outfile = open('query_filtered_clustering_blast_results', 'w')
twoway_filter_outfile =  open('twoway_filtered_clustering_blast_results', 'w')
for i, line in enumerate(blast_results):
    linespl = line.strip().split()
    if float(linespl[2]) <= 30.:
        low_identity_discarded += 1
        continue
    qlen = query_length[linespl[0]]
    if float(linespl[3])/qlen >= 0.5:
        query_filter_outfile.write(line)
        slen = query_length[linespl[1]]
        if float(linespl[3])/slen >= 0.5:
            twoway_filter_outfile.write(line)
        else:
            twoway_filter_discarded += 1
    else:
        query_filter_discarded += 1
        twoway_filter_discarded += 1
    if not i % 1e08:
        print('Processed %i lines' % i)
        print('%i low identity HSPs discarded' % low_identity_discarded)
        print('%i low query cover discarded' % query_filter_discarded)
        print('%i low twoway cover discarded' % twoway_filter_discarded)


#print('Low-identity HSPs removed; %i lines remaining' % len(blast_results))
#query_filtered = [l for l in blast_results if int(l[3])/query_length[l[0]] >= 0.5]
#print('Query filter completed; %i lines remaining' % len(query_filtered))
#twoway_filtered = [l for l in query_filtered if int(l[3])/query_length[l[1]] >= 0.5]
#print('Two-way filter completed; %i lines remaining' % len(twoway_filtered))
#print('Saving query filter results')
#query_filtered = '\n'.join(['\t'.join(l) for l in query_filtered]) + '\n'
                

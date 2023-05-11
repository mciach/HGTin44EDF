# HGTin44EDF
Horizontal Gene Transfer in 44 Early Diverging Fungi

This repository contains a detailed description of the pipeline from (Ciach et al., 2023) to detect horizontally transfered protein-encoding genes in a set of proteomes.    
The individual steps are described below. The repository also contains necessary Jupyter Notebooks and Python3 scripts to perform the analysis.  

 ![The flow chart of the pipeline](Figures/flowchart_nowy.png)

# Required data

1. Directory `Proteomes/` containing FASTA files with target proteomes. 
   FASTA headers should ideally contain only the accession of the protein.  
2. Directory `GFFs/` containing GFF files corresponding to the proteomes.   
3. Additional files (supplied with the protocol):
    1. `org2taxid.tsv`: A tab-separated file containing proteome filenames, organism names, species TaxIDs and family TaxIDs. Supplied with the repository in the Additional_data directory. For other studies, this needs to be prepared manually.   
4. Additional files (require manual preparation):
    1. `acc2taxid.txt`:  A tab-separated file containing all protein accessions from `Proteomes/` and the corresponding species TaxID. Easy to generate using `org2taxid.tsv`. 
    2. `acc2filename.txt`: A tab-separated file containing all protein accessions from `Proteomes/` and the corresponding proteome filename. Easy to generate using `org2taxid.tsv`. 
    3. `basal_accessions.txt`: A list of all target accessions (essentially, the first column of `acc2taxid.txt`).  
    4. `all_proteins.fa`: a concatenated FASTA of all proteomes. Can be generated with `cat Proteomes/* > all_proteins.fa`.  

# Data preparation (workflow step 1):

1. Generate a tab-separated file `acc2scaffold.tsv` mapping the acessions of protein encoding genes to their genes' scaffolds.  
BASH commands:   
```
cat GFFs/*.gff | cut -f1,9 | grep 'CDS' | grep 'protein_id' | sed 's/\t.*=/\t/' > tmp_acc2scaffold.tsv
sort -uo tmp_acc2scaffold.tsv tmp_acc2scaffold.tsv  
paste <(cut -f2 tmp_acc2scaffold.tsv) <(cut -f1 tmp_acc2scaffold.tsv) > acc2scaffold.tsv
rm tmp_acc2scaffold.tsv
```

2. Mask low-complexity regions.  
BASH command using ncbi-seg:
```
ncbi-seg all_proteins.fa 12 1.8 2.0 -x X > masked_proteins.fa
```
*Note: SEG parameters taken from http://manpages.ubuntu.com/manpages/bionic/man1/ncbi-seg.1.html, 
giving less masking, so that we'll discard less proteins. This is optional, default parameters can also be used.*

3. Remove heavily masked proteins, short/long proteins, remove atypical characters. 
BASH command using a supplied Python script `standardize_proteins.py`:
```
python3 standardize_proteins.py masked_proteins.fa standardized_proteins.fa
```
*Note: The resulting file `standardized_proteins.fa` is now the main FASTA file that will be screened for HGT.*

# Data processing 
## 1. Discarding ancestrally fungal proteins and major contaminants (workflow steps 2 and 3):
1. Create a custom BLASTp database `target_fungi` with target proteins.  
BASH command using ncbi blast+:
```
    makeblastdb -in ../standardized_proteins.fa -dbtype prot -out target_fungi -parse_seqids -taxid_map ../acc2taxid.txt
```

2. All-vs-all BLASTp to asses the taxonomic distribution of the target proteins across the target fungi.  
BASH command using ncbi blast+:
```
blastp -db target_fungi -query standardized_proteins.fa -out initial_blast_results -evalue 0.00001 -outfmt "6 qaccver saccver staxid" -num_threads 60 -qcov_hsp_perc 50
```
*Note: By their definition, E-values are proportional to the target data base size. 
As a consequence, BLASTing against smaller data bases always yields smaller E-values than against large data bases. 
To take this fact into account, in this step we make a rough adjustment of our e-value threshold and set it to 0.00001. 
This is mostly optional - if you use a higher e-value, you'll increase the sensitivity, but decrease the specificity of this BLASTp search.*

3. Parse the BLASTp results.  
Use the supplied Jupyter Notebook: "1. Discarding ancestral sequences and major contaminants.ipynb".   
Implemented steps:  
    1. Select proteins with homology confined to a single fungal family (using the org2taxid.tsv table). 
    2. Discard proteins on contigs without any homology to any other fungal proteome.  
    3. Save the selected proteins in FASTA format in the file `global_blast_query.fa.`


## 2. Performing a BLASTp of the selected target proteins against the NR data base and filtering the results ("global BLASTp"; workflow step 4, 5)
1. Run a BLASTp of the global_blast_query.fa against a local copy of the NR data base.   
Set a query cover filter so that we don't get spurious homologies caused by small highly conserved regions (this is very prevalent in these kinds of studies, because HGT is often associated with gene fusion).  
Example command (it is recommended, however, to parallelize the computations, see e.g. [this link](https://bioinformaticsworkbook.org/dataAnalysis/blast/running-blast-jobs-in-parallel.html#gsc.tab=0)):
```
blastp -db path_to_nr -query global_blast_query.fa -out global_blast_results -evalue 0.001 -outfmt "6 qaccver saccver staxid" -num_threads 60 -qcov_hsp_perc 50
```

2. Extract the accessions of all homologs and create a homolog accession - taxid mapping.   
Remove sequences with taxid=0 (unknown taxon). 
BASH commands:
```
cut -f2,3 global_blast_results | sort -u > homolog_accession_taxid_table
grep -Pv '\t0$' homolog_accession_taxid_table > filtered_homolog_taxid_table
```

3. Filter the results - discard viral and environmental sequences 
BASH commands using the supplied Python script `remove_viral_homologs.py`:   
```
python3 remove_viral_homologs.py filtered_homolog_taxid_table final_homolog_taxid_table
cut -f1 final_homolog_taxid_table > final_accessions
```

4. Fetch the sequences of homologs.   
BASH commands:
```
blastdbcmd -db path_to_nr -dbtype prot -entry_batch final_accessions -target_only > final_homolog_sequences.fa  
```
*Note: -target_only is crucial here, otherwise you'll get a few hundred synonyms per sequence.* 

5. Mask the low-complexity regions of homologs.
BASH command:  
```
ncbi-seg -in final_homolog_sequences.fa -infmt fasta -outfmt fasta -x X -out masked_final_sequences.fa
```

6. Discard too short, too long, or too masked sequences.
BASH command using a supplied Python script `standardize_proteins.py`:
```
python3 standardize_proteins.py masked_final_sequences.fa standardized_final_sequences.fa
```
*Note: `standardized_final_sequences.fa` is now the main file that will be used for subsequent computations.*  
*Note 2: You may want to adjust the script to relax the masking threshold - e.g. allow for up to 70% masking, to retain more homologs.  
For subsequent computations, this threshold is less important than the length threshold, and may be more freely adjusted according to your needs.  
For example, if you encounter preoblems with alignment quality or tree inference, it may be worth to make this threshold more stringent to remove low-complexity proteins, which are difficult to align.  
If, on the other hand, you encounter problems with clustering because your clusters are too small, and you need to retain more sequences, it may be worth to relax the masking threshold. It's a trade-off, as always.  
In our case, a 70% masking threshold (i.e. setting `low_complexity = aa_nbs/lengths <= 0.3` in line 59 of the script) worked well.*

7. Generate a tsv file with protein sequence lengths called `sequence_length_table.txt`. The file needs to have two tab-separated columns, with protein accession in the first one and protein sequence length in second one. Easy to do manually e.g. using `Biopython`, so no script or command is supplied here.  

## 3. Clustering sequences - fist stage (workflow step 6)
1. Create a BLASTp database from `standardized_final_sequences.fa` for all-vs-all BLASTp.  
BASH commands:  
```
makeblastdb -in standardized_final_sequences.fa -dbtype prot -out clustering_dbase
```

2. Run BLASTp.  
BASH commands:  
```
blastp -db clustering_dbase -out clustering_blast_results -word_size 6 -threshold 21 -evalue 1e-06  -outfmt \"6 qaccver saccver pident length evalue bitscore\" -query standardized_final_sequences.fa
```
*Note 1: Here we adjust the word size and score threshold to speed up the computation. This makes BLASTp less sensitive, but more specific, which we can afford in this analysis as we're not looking for distant homologs.   
Note 2: We make a rough adjustment of the e-value threshold for the same reason as in the first BLASTp.*

3. Filter BLASTp results: select only homolog pairs with identity >= 30%, query cover >= 50%, and subject cover >= 50%.   
BLASTp HSPs not meeting these criteria are discarded (i.e. we assume there is no homology).  
This makes the homology graph for MCL more sparse, and leads to more clusters which are also more biologically plausible.   
Otherwise, clusters are joined by small, conserved regions with low e-values.     
Effectively, this step increases the specificity of homology search, but decreases its sensitivity.  
BASH command using a supplied Python script `filter_blast_results.py`:  
```python filter_blast_results.py```
*Note: `filter_blast_results.py` reads results from clustering_blast_results and generates two files: `query_filtered_clustering_blast_results` and `twoway_filtered_clustering_blast_results`, for two filtering strategies. The file `twoway_filtered_clustering_blast_results` is the one with full filtering, used in subsequent steps; `query_filtered_clustering_blast_results` can be discarded or used to compare the results.*   
*Note 2: `filter_blast_results.py` uses data from the `sequence_length_table.txt` file created in step 2.7.*

4. Cluster the file `twoway_filtered_clustering_blast_results` with MCL, inflation parameter = 1.7.   
Instructions available at https://micans.org/mcl/.  
Save the clustering results in a file `out.mcl_twoway_filter.mci.I17`.  
*Note: The precise value of the inflation parameter is of a minor importance; For example, the sometimes recommended value of 1.8 gives similar results.   
There is a specificity/sensitivity trade-off here just as in the previous steps, which does not have a major influence on the final biological conclusions, but may result in different 
numbers of detected HGTs and different numbers of false positive results. It is recommended to inspect the results after this step, e.g. by selecting random clusters and comparing with online BLASTp search to see if the homologous proteins are within the sequence cluster, and whether there are spurious sequences in the cluster.*  


## 4. Cluster processing and filtering (workflow step 7)
1. Use the Jupyter Notebook 2. "Selecting raw clusters (all stages)" to process clusters from the file `out.mcl_twoway_filter.mci.I17`.  
The notebook will save the FASTA files with sequences from each cluster in the directory `first_round_raw_clusters`.  
Note: The notebook can be used to filter clusters in this step as well as after re-clustering in workflow step 8 by modifying the value of the variable `cluster_file` in the notebook.  
Applied filters:
    1. Discard clusters without any target protein.
    2. Discard clusters with fungal proteins from more than one taxonomic fungal phylum.
   This step is optional. It discards clusters which are spread throughout fungi, in order to limit the risk of false positives
   and to speed up the computations. It reduces the sensitivity of the analysis (we'll obtain less HGTs in the end), 
   but also removes clusters which may be difficult to analyze and therefore may introudce errors.  
    3. Discard clusters where non-fungal species constitute less than 60% of all species.
   This step is optional - it simply discards clusters that are unlikely to contain HGTs into fungi
   in order to speed up computations. These clusters are unlikely to contain HGTs into the fungi
   because the taxonomic composition is mostly fungal.



## 5. Second-stage contaminant filtering (workflow step 8)
1. Use the supplied Jupyter Notebook 3. "Contaminant screening and removal", section "Generating contaminant filtering blast query",  
to select a random sample of proteins for a BLASTp agains the NR data base.  
For each target protein, the notebook selects up to 10 proteins encoded on the same contig.   
The notebook saves the selected BLASTp query in a directory called `contaminant_filtering_blast`, in the subdirectory `contaminant_filtering_blast/sequences`.   
*Note 1: In principle, this step can be done at any earlier stage, even at the beginning of the analysis.   
However, it would be very computationally costly. Performing it at this stage, after an initial filtering of clusters, saves a lot of computational time.   
Note 2: Selecting up to 10 proteins is quite arbitary, the more the better (but again, costlier).*    

2. Run a BLASTp of the contaminant filtering query against the NR data base.  
Example BASH command with parallelized BLASTp (evoked in the `contaminant_filtering_blast` directory):  
```
mkdir blast 
parallel --jobs=20 blastp -db path_to_local_NR -query sequences/{1} -out blast/{1}.blast -outfmt \"6 qseqid sseqid pident length evalue staxid\" -num_threads 4 -word_size 6 -threshold 21 -evalue 1e-06 -max_target_seqs 100 ::: $(ls sequences)   
cat blast/* > blast_results
```
*Note: the local NR data base needs to be constructed with TaxID information in order for the `staxid` keyword to work.*

3. Use the supplied Jupyter Notebook 3. "Contaminant screening and removal", section "Detecting contaminants", 
to parse the results of BLASTp from the previous step and to discard the target proteins for which we didn't detect any protein with a fungal first hit agains NR
encoded in the same contig. These proteins are removed from clusters.  
Then, the notebook repeats filtering of clusters as in the previous step of the workflow.  
This step is optional, but at this point we may get some clusters without target proteins, so their further processing is unneccessary and they can be safely discarded.  
The notebook will then generate filtered cluster FASTAs and save them in a directory called `first_round_filtered_clusters`.     


## TBC
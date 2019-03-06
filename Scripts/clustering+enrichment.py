# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 14:02:38 2019

@author: Milica

This script performs afiinity propagation and an enrichment analysis.
"""

import os
import numpy as np
import pandas as pd
from sklearn.cluster import AffinityPropagation
from orangecontrib.bioinformatics import go
import mygene

os.chdir('F:\\PDS\\Project')

sim_mat = np.loadtxt('Data\\adjecency_matrix.txt', delimiter=',') # + and - regulation 
                                                                  
af = AffinityPropagation(affinity='precomputed').fit(sim_mat)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_
n_clusters_ = len(cluster_centers_indices)

# Average cluster size

lab = np.repeat(0, n_clusters_)
for i in labels:
    lab[i] += 1

np.mean(lab)
np.median(lab)

# _____________________________________________________________________________

# Makes a list of all genes and converts gene symbols to Entrez IDs

data = pd.read_csv('sgd.txt', sep='\t', na_values = '?', dtype='float')
gene_labels = list(data.columns)

mg = mygene.MyGeneInfo()
results = mg.querymany(gene_labels, scopes='symbol', species=559292, as_dataframe=True)

# _____________________________________________________________________________

def get_current_cluster(positions, all_genes):
    # takes genes form all_genes from positions and makes a list
    gene_list = []
    for value in positions:
        gene_list.append(all_genes[int(value)])
    return gene_list

# Enrichment analysis 
    
ontology = go.Ontology()
annotations = go.Annotations("4932", ontology=ontology)
annotations.ontology.set_slims_subset('goslim_yeast')
enriched = pd.DataFrame(columns=['GenesInCluster', 'EntrezID', 'EnrichedGO'], 
                        index=range(n_clusters_))

for i in range(n_clusters_):
    positions = np.argwhere(labels==i) #returns indexes of genes in cluster i
    current_cluster = get_current_cluster(positions, gene_labels) 
    enriched.GenesInCluster[i] = current_cluster
    enriched.EntrezID[i] = mg.querymany(current_cluster, 
                scopes='symbol', species=559292, as_dataframe=True).entrezgene
    enriched.EnrichedGO[i] = annotations.get_enriched_terms(enriched.EntrezID[i], 
                slims_only=True)

enriched.to_csv('Data\\enriched.csv', sep=',', header=True)

# _____________________________________________________________________________
    
hits = pd.DataFrame(columns=('ClusterGenes', 'HitGenes', 'GO_term', 'depth', 'p_value'), 
                    index=range(327))

i = 0
for enr_dict in enriched.EnrichedGO:
    arg = int(np.where(enriched.EnrichedGO == enr_dict)[0])
    for go_id, (genes, p_value, ref) in enr_dict.items():
        if p_value < 0.1:
            i += 1
            hits.ClusterGenes[i] = enriched.GenesInCluster[arg]
            hits.HitGenes[i] = genes
            hits.GO_term[i] = go_id
            hits.depth[i] = ontology.term_depth(go_id)
            hits.p_value[i] = p_value

len(np.unique(hits['ClusterGenes']))

hits.to_csv('Data\\hits.csv', sep=',', header=True)    
           
# _____________________________________________________________________________

hits = pd.read_csv('Data\\hits.csv', sep=',', index_col=0)

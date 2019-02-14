# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 19:44:15 2019

@author: Milica
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from scipy import stats
from matplotlib import pyplot as plt

os.chdir('F:\\PDS\\Project')

# Read the data
data = pd.read_csv('sgd.txt', sep='\t', na_values = '?', dtype='float')
labels = list(data.columns)

# Replace missing values with gene median (column)??
for i in range(data.shape[1]):
    data.iloc[:, i] = data.iloc[:, i].fillna(np.nanmedian(data.iloc[:, i]))

# Store values as numpy matrix
gene_expr = np.matrix(data.values)

# Compute correlation matrix
corr_mat = np.zeros((1799,1799))
for i in range(1799):
    for j in range(1799):
        corr_mat[i, j] = stats.pearsonr(gene_expr[:,i], gene_expr[:,j])[0]

# Export original correlation matrix (before the permutation)
corr_org = corr_mat.reshape(1799*1799,)
np.savetxt('original_corr_mat.txt', corr_mat, delimiter=',') 

# Present original matrix as heatmap
fig1 = sns.heatmap(corr_mat)
fig1 = fig1.get_figure()
fig1.savefig("Figures//original_corr_mat_HM.png", dpi = 900)

# Remove the noise

perm_mat = gene_expr
corrs = np.array(np.repeat(0, 10000), dtype = 'float')

for i in range(10000):    
    for j in range(740):
        perm_mat[j,:] = np.random.permutation(perm_mat[j,:])    
    a = perm_mat[:,np.random.randint(1799, size=1)]
    b = perm_mat[:,np.random.randint(1799, size=1)]
    corrs[i] = stats.pearsonr(a, b)[0]
   
np.savetxt('permuted_corr_mat.txt', corrs, delimiter=',')

# Present correlations with histograms, treba dodati legende

fig2 = plt.figure(figsize=(8, 16), dpi = 150)
fig2.add_subplot(211)
sns.distplot(corrs, hist=True, kde=True, rug=False, color='blue')
fig2.add_subplot(212)
sns.distplot(corr_org, hist=True, kde=True, rug=False, color='green')
plt.show()
fig2.savefig('Figures\\correlation_histograms.png', dpi = 150)
plt.clf()

# Set treshold
#only gives positive correlation? ? !!!
treshold_up = np.percentile(corrs, 99)
treshold_down = np.percentile(corrs, 1)


# Make weighted graph
G = nx.Graph()

for i in range(1799):
    for j in range(i + 1, 1799):
        if corr_mat[i, j]>treshold_up or corr_mat[i, j]<treshold_down:
            G.add_edge(labels[i], labels[j], weight = np.abs(corr_mat[i, j]))

# Export as weighted edgelist
nx.readwrite.edgelist.write_weighted_edgelist(G, 'G.edgelist')

# Degree distribution
sns.distplot(list(d for n, d in G.degree()), hist=True, kde=False, rug=False)


# Division into subgraphs

global_props = pd.DataFrame(np.zeros((17, 6)), dtype='float', 
                            columns=['Order', 'Size', 'Radius', 
                                     'Diameter', 'Clust_coeff', 'Avg_path_len'])
cnt = 0
for subgraph in [G.subgraph(c) for c in nx.connected_components(G)]:
    nx.readwrite.edgelist.write_weighted_edgelist(subgraph, 'Graphs//G' + str(cnt) + '.edgelist')
    global_props['Order'][cnt] = subgraph.order()
    global_props['Size'][cnt] = subgraph.size()
    global_props['Radius'][cnt] = nx.radius(subgraph)
    global_props['Diameter'][cnt] = nx.diameter(subgraph)
    global_props['Clust_coeff'][cnt] = nx.average_clustering(subgraph)
    global_props['Avg_path_len'][cnt] = nx.average_shortest_path_length(subgraph)
    cnt += 1

global_props.to_csv('subgraphs_properties.csv')

"""
# Old Subgraphs with more that 2 nodes (nodes/edges)
G0 = G.subgraph(list(nx.connected_components(G))[0])  # 11/18   
G1 = G.subgraph(list(nx.connected_components(G))[1])  # 830/14378    
G5 = G.subgraph(list(nx.connected_components(G))[5])  # 5/5
G7 = G.subgraph(list(nx.connected_components(G))[7])  # 8/13  
"""
G0 = G.subgraph(list(nx.connected_components(G))[0])  # 903/32482   
G5 = G.subgraph(list(nx.connected_components(G))[5])  # 5/5    
G6 = G.subgraph(list(nx.connected_components(G))[6])  # 3/2

nx.readwrite.edgelist.write_weighted_edgelist(G0, 'Graphs_new\\G0.edgelist')
nx.readwrite.edgelist.write_weighted_edgelist(G5, 'Graphs_new\\G5.edgelist')
nx.readwrite.edgelist.write_weighted_edgelist(G6, 'Graphs_new\\G6.edgelist')

"""
# G1 the bggest graph - other properties and node ranking

# Degree distribution
sns.distplot(list(d for n, d in G1.degree()), hist=True, kde=False, rug=False)

# Shortest paths distribution
G1_shorthest_paths = list(dict(nx.shortest_path_length(G1)).values())
G1_paths = []
for node_paths in G1_shorthest_paths:
    G1_paths += list(node_paths.values())

sns.distplot(G1_paths, hist=True, kde=False, rug=False)

# Plot G1
fig3 = plt.figure(figsize = (20, 20))
pos = nx.fruchterman_reingold_layout(G1)
nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color='b')
plt.title('Subgraph G1')
plt.show()
fig3.savefig('Figures//G1.png')

# Centrality measures
G1_centrality = pd.DataFrame(np.zeros((830, 5)), dtype='float', 
                            columns=['Degree', 'Eccentricity', 'Closness', 
                                     'Betweeness', 'PageRank'])

G1_centrality['Degree'] = [value for value in nx.degree_centrality(G1).values()]
G1_centrality['Eccentricity'] = [1/value for value in nx.eccentricity(G1).values()]
G1_centrality['Closness'] = [value for value in nx.closeness_centrality(G1).values()]
G1_centrality['Betweenes'] = [value for value in nx.betweenness_centrality(G1).values()]
G1_centrality['PageRank'] = [value for value in nx.pagerank(G1).values()]

G1_max = {}
for column in G1_centrality:
    G1_max[column] = np.argsort(G1_centrality[column])[0:10]

fig4 = plt.figure(figsize = (20, 20))
pos = nx.fruchterman_reingold_layout(G1)
#node_sizes = [d*10 for d in degrees]
node_color = np.repeat(1, 830)
node_color[list(G1_max['PageRank'])] = 5
nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
plt.show()
fig4.savefig('G1_PageRank.png')


fig4 = plt.figure(figsize = (20, 20))
pos = nx.kamada_kawai_layout(G1)
#node_sizes = [d*10 for d in degrees]
node_color = np.repeat(1, 830)
node_color[list(G1_max['PageRank'])] = 5
nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
plt.show()
fig4.savefig('G1_PageRank_KK.png')


for key in G1_max:
    fig = plt.figure(figsize = (20, 20))
    pos = nx.fruchterman_reingold_layout(G1)
    node_color = np.repeat(1, 830)
    node_color[list(G1_max[key])] = 5
    nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
    nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
    plt.show()
    fig.savefig('G1_' + key + '.png')

"""
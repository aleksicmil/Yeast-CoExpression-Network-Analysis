# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 18:58:35 2019

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

G1 = nx.nx.readwrite.edgelist.read_weighted_edgelist('Graphs\\G1.edgelist')

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

# Dictonary of 10 nodes with highest centrality measures
G1_max = {}
for column in G1_centrality:
    G1_max[column] = np.argsort(G1_centrality[column])[0:10]


# Plot highest page rank
fig4 = plt.figure(figsize = (20, 20))
pos = nx.fruchterman_reingold_layout(G1)
#node_sizes = [d*10 for d in degrees]
node_color = np.repeat(1, 830)
node_color[list(G1_max['PageRank'])] = 5
nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
plt.show()
fig4.savefig('G1_PageRank.png')

# Plot highest page rank
fig4 = plt.figure(figsize = (20, 20))
pos = nx.kamada_kawai_layout(G1)
#node_sizes = [d*10 for d in degrees]
node_color = np.repeat(1, 830)
node_color[list(G1_max['PageRank'])] = 5
nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
plt.show()
fig4.savefig('G1_PageRank_KK.png')

# Plots G1 wih G1_max nodes in different color
for key in G1_max:
    fig = plt.figure(figsize = (20, 20))
    pos = nx.fruchterman_reingold_layout(G1)
    node_color = np.repeat(1, 830)
    node_color[list(G1_max[key])] = 5
    nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
    nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
    plt.show()
    fig.savefig('G1_' + key + '.png')


def get_genes(graph, max_values):
    node_names = [node for node in graph.nodes]
    genes = []
    for rank in max_values:
        genes.append(node_names[rank])
    return genes

        
betweenes_genes = get_genes(G1, G1_max['Betweeness'])
file = open('Data\\Betweeness_genes.tab', 'w')
for gene in betweenes_genes: file.write(gene + '\n')
file.close()


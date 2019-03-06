# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:58:25 2019

@author: Milica

G is the biggest subgraph from the network made of *only* positive gene regulation
relationships.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
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
fig3.savefig('Plots//G1.png')

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
node_color = np.repeat(1, 830)
node_color[list(G1_max['PageRank'])] = 5
nx.draw_networkx_edges(G1, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G1, pos, alpha=0.6, node_color=node_color)
plt.show()
fig4.savefig('G1_PageRank.png')


fig4 = plt.figure(figsize = (20, 20))
pos = nx.kamada_kawai_layout(G1)
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


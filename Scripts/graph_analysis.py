# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 22:57:08 2019

@author: Milica

This script takes adjecency matrix and makes a weighted graph. The graph is further 
separated into connected components.

The largest component (G0) is analysed and plotted. 

'Data\\adjecency_matrix.txt' includes *both* positive and negative gene regulation
relationships.

"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from matplotlib import pyplot as plt

os.chdir('F:\\PDS\\Project')

# Read adjecency matrix produced after cleaning initial data
adj_mat = np.loadtxt('Data\\adjecency_matrix.txt', delimiter=',')

# Make weighted graph
G = nx.from_numpy_matrix(adj_mat)
nx.readwrite.edgelist.write_weighted_edgelist(G, 'Graphs_new\\G.edgelist')

## _________________
## Graph description
## _________________

# Degree distirubtion
f_deg = sns.distplot(list(d for n, d in G.degree()), hist=True, kde=False, rug=False)
f_deg = f_deg.get_figure()
f_deg.savefig('Plots//G_degree_distribution.png', dpi = 100)

# Shortest paths distribution
G_shorthest_paths = list(dict(nx.shortest_path_length(G)).values())
G_paths = []
for node_paths in G_shorthest_paths:
    G_paths += list(node_paths.values())

f_paths = sns.distplot(G_paths, hist=True, kde=False, rug=False).get_figure()
f_paths.savefig('Plots\\G_shortest_path_distribtion.png', dpi = 100)


## -----------------------
## Division into subgroups
## -----------------------


# Produces and exports a dataframe with properties of different subgraphs 

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

global_props.to_csv('Data\\subgraphs_properties.csv')

#global_props = pd.read_csv('subgraphs_properties.csv', sep=',', index_col=0)
# Biggest subgraphs (number of nodes/number of enges)

G0 = G.subgraph(list(nx.connected_components(G))[0])  # 903/32482   
G5 = G.subgraph(list(nx.connected_components(G))[5])  # 5/5    
G6 = G.subgraph(list(nx.connected_components(G))[6])  # 3/2

# Exports biggest graphs as edgelists

nx.readwrite.edgelist.write_weighted_edgelist(G0, 'Graphs_new\\G0.edgelist')
nx.readwrite.edgelist.write_weighted_edgelist(G5, 'Graphs_new\\G5.edgelist')
nx.readwrite.edgelist.write_weighted_edgelist(G6, 'Graphs_new\\G6.edgelist')

## ---------------------------------------------
## G0 - biggest subgraph from pos+neg adj matrix
## ---------------------------------------------


G0 = nx.nx.readwrite.edgelist.read_weighted_edgelist('Graphs_new\\G0.edgelist')


# Degree distribution

sns.distplot(list(d for n, d in G0.degree()), hist=True, kde=False, rug=False)

# Shortest paths distribution

G0_shorthest_paths = list(dict(nx.shortest_path_length(G0)).values())
G0_paths = []
for node_paths in G0_shorthest_paths:
    G0_paths += list(node_paths.values())

sns.distplot(G0_paths, hist=True, kde=False, rug=False)

# Plot G0
fig3 = plt.figure(figsize = (20, 20))
pos = nx.fruchterman_reingold_layout(G0)
nx.draw_networkx_edges(G0, pos, alpha=0.2, edge_color='k')
nx.draw_networkx_nodes(G0, pos, alpha=0.6, node_color='b')
plt.title('Subgraph G0')
plt.show()
fig3.savefig('Figures//G0.png')

# Centrality measures

G0_centrality = pd.DataFrame(np.zeros((903, 5)), dtype='float', 
                            columns=['Degree', 'Eccentricity', 'Closness', 
                                     'Betweeness', 'PageRank'])

G0_centrality['Degree'] = [value for value in nx.degree_centrality(G0).values()]
G0_centrality['Eccentricity'] = [1/value for value in nx.eccentricity(G0).values()]
G0_centrality['Closness'] = [value for value in nx.closeness_centrality(G0).values()]
G0_centrality['Betweenes'] = [value for value in nx.betweenness_centrality(G0).values()]
G0_centrality['PageRank'] = [value for value in nx.pagerank(G0).values()]
G0_centrality.index = nx.degree_centrality(G0).keys()

# Dictonary of 10 nodes with highest centrality measures

G0_max = {}
for column in G0_centrality:
    size = len(G0_centrality[column])
    current_values = np.array(G0_centrality[column])
    G0_max[column] = G0_centrality.index[np.argsort(current_values)[(size-11):(size-1)]]

gene_names = {}
for i in range(len(G0.nodes)):
    gene_names[np.array(G0.nodes)[i]] = i
    

# Plot highest page rank, fruchtman reingold layout
    
fig4 = plt.figure(figsize = (20, 20))
pos = nx.fruchterman_reingold_layout(G0)
node_color = np.repeat('k', 903)
for gene in G0_max['PageRank']:
    node_color[gene_names[gene]] = 'b'
nx.draw_networkx_edges(G0, pos, alpha=0.1, edge_color='k')
nx.draw_networkx_nodes(G0, pos, alpha=0.6, node_color=node_color, node_size=100)
plt.show()
fig4.savefig('Plots\\G0_PageRank_new.png')



# Plots G0 wih G0_max nodes in different color

for key in G0_max:
    fig = plt.figure(figsize = (20, 20))
    pos = nx.fruchterman_reingold_layout(G0)
    node_color = np.repeat('k', 903)
    for gene in G0_max[key]:
        node_color[gene_names[gene]] = 'b'
    nx.draw_networkx_edges(G0, pos, alpha=0.1, edge_color='k')
    nx.draw_networkx_nodes(G0, pos, alpha=0.6, node_color=node_color,node_size=90)
    plt.show()
    fig.savefig('Plots\\G0_' + key + '.png')

# Ego graphs

ego = nx.ego_graph(G0, G0_max['Betweeness'][9], radius=1)
fig = plt.figure(figsize = (10, 10))
pos = nx.fruchterman_reingold_layout(ego)
nx.draw_networkx_edges(ego, pos, alpha=0.1, edge_color='k')
nx.draw_networkx_nodes(ego, pos, alpha=0.6, node_size=90, node_color='b')
nx.draw_networkx_labels(ego, pos, labelsls=ego.nodes)
fig.savefig('Plots\\Ego' + 'Betweenes9' + '.png')



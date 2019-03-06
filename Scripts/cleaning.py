# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 19:44:15 2019

@author: Milica

Script that reads in the data, fills in the missing values, removes the noise, 
and makes the adjecency matrix for the future analysis.

"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from matplotlib import pyplot as plt

os.chdir('F:\\PDS\\Project')

# Read the data

data = pd.read_csv('Data\\sgd.tab', sep='\t', na_values = '?', dtype='float')
labels = list(data.columns)


# Replace missing values with gene(column) median and store it as numpy matrix

for i in range(data.shape[1]):
    data.iloc[:, i] = data.iloc[:, i].fillna(np.nanmedian(data.iloc[:, i]))

gene_expr = np.matrix(data.values)

# Compute correlation matrix and export it

corr_mat = np.zeros((1799,1799))
for i in range(1799):
    for j in range(1799):
        corr_mat[i, j] = stats.pearsonr(gene_expr[:,i], gene_expr[:,j])[0]

corr_org = corr_mat.reshape(1799*1799,)
np.savetxt('Data\\original_corr_mat.txt', corr_mat, delimiter=',') 

# Present original matrix as heatmap

fig1 = sns.heatmap(corr_mat)
fig1 = fig1.get_figure()
fig1.savefig("Plots//original_corr_mat_HM.png", dpi = 900)


# Remove the noise
perm_mat = gene_expr
corrs = np.array(np.repeat(0, 10000), dtype = 'float')

for i in range(10000):    
    for j in range(740):
        perm_mat[j,:] = np.random.permutation(perm_mat[j,:])    
    a = perm_mat[:,np.random.randint(1799, size=1)]
    b = perm_mat[:,np.random.randint(1799, size=1)]
    corrs[i] = stats.pearsonr(a, b)[0]
   
np.savetxt('Data\\permuted_corr_mat.txt', corrs, delimiter=',')

# Present correlations with histograms, treba dodati legende

fig2 = plt.figure(figsize=(8, 16), dpi = 150)
fig2.add_subplot(211)
sns.distplot(corrs, hist=True, kde=True, rug=False, color='blue')
fig2.add_subplot(212)
sns.distplot(corr_org, hist=True, kde=True, rug=False, color='green')
plt.show()
fig2.savefig('Plots\\correlation_histograms.png', dpi = 150)
plt.clf()

# Set treshold and create an adjecency matrix

corr_mat = np.loadtxt('Data\\original_corr_mat.txt', delimiter=',')
corrs = np.loadtxt('Data\\permuted_corr_mat.txt', delimiter=',')

treshold_up = np.percentile(corrs, 99)
treshold_down = np.percentile(corrs, 1)


adj_mat = np.copy(corr_mat) # both positive and negative correlations
for i in range(1799):
    for j in range(1799):
        if np.logical_or(adj_mat[i, j] > treshold_up, adj_mat[i, j] < treshold_down):
            adj_mat[i, j] = np.abs(adj_mat[i, j])
        else:
            adj_mat[i, j] = 0


neg_adj_mat = np.copy(corr_mat) # only negative correlations
for i in range(1799):
    for j in range(1799):
        if neg_adj_mat[i, j] < treshold_down:
            neg_adj_mat[i, j] = np.abs(adj_mat[i, j])
        else:
            neg_adj_mat[i, j] = 0

np.fill_diagonal(neg_adj_mat, 1)

pos_adj_mat = np.copy(corr_mat) # only positive correlations
for i in range(1799):
    for j in range(1799):
        if pos_adj_mat[i, j] > treshold_up:
            pos_adj_mat[i, j] = pos_adj_mat[i, j]
        else:
            pos_adj_mat[i, j] = 0


np.savetxt('Data\\adjecency_matrix.txt', adj_mat, delimiter=',')
np.savetxt('Data\\negative_adjecency_matrix.txt', neg_adj_mat, delimiter=',')
np.savetxt('Data\\positive_adjecency_matrix.txt', pos_adj_mat, delimiter=',')


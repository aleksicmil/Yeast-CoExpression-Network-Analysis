# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 14:02:38 2019

@author: Milica
"""

import os
import numpy as np
#import pandas as pd
#import seaborn as sns
#import networkx as nx
#from matplotlib import pyplot as plt
from sklearn.cluster import AffinityPropagation
from sklearn import metrics

os.chdir('F:\\PDS\\Project')

sim_mat = np.loadtxt('Data\\adjecency_matrix.txt', delimiter=',')

af = AffinityPropagation(preference=-50).fit(X)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_



# Yeast-CoExpression-Network-Analysis
## Project for course "Programming for Data Science"

This is a project done for the course Programming for Data Science, by professor Sanja Brdar, 
on the master program "Data Sciencce" at University of Novi Sad, Faculty of Sciences. The project task is 
described in the [Project 6 - Bioinformatics Data Analysis.pdf](https://github.com/aleksicmil/Yeast-CoExpression-Network-Analysis/blob/master/Project%206%20-%20Bioinformatics%20Data%20Analysis.pdf), 
and in the [project presentation](https://github.com/aleksicmil/Yeast-CoExpression-Network-Analysis/blob/master/Presentation/Aleksic_PDS_yeast.pptx) you 
can find more detailed explanation of my work.  

Work is divided into four scripts:
  * [cleaning.py](https://github.com/aleksicmil/Yeast-CoExpression-Network-Analysis/blob/master/Scripts/cleaning.py) - Script that 
reads in the data, fills in the missing values, removes the noise, 
and makes the adjecency matrix for the future analysis
  * [graph_analysis.py](https://github.com/aleksicmil/Yeast-CoExpression-Network-Analysis/blob/master/Scripts/graph_analysis.py) - This 
script takes adjecency matrix and makes a weighted graph. The graph is further 
separated into connected components. The largest component is further analysed. Both positive and negative correlations between
genes are included.
  * [G1_analysis.py](https://github.com/aleksicmil/Yeast-CoExpression-Network-Analysis/blob/master/Scripts/G1_analysis.py)- G1 is 
the biggest subgraph from the network made of __only__ positive gene regulation
relationships.
  * [clustering+enrichment.py](https://github.com/aleksicmil/Yeast-CoExpression-Network-Analysis/blob/master/Scripts/clustering+enrichment.py) - This 
script clusters genes using affinity propagation and performs an enrichment analysis.
 
**Data** is from [yeastgenome.org](www.yeastgenome.org) and contains information on 1799 genes and their expression in 740 different experimental contitions.
More information about this can be found [here](https://sites.google.com/view/yeastgenome-help/function-help/expression-data).

Plots and constructed graphs are grouped in the folders with the same name. 
 

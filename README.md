# disease-drug-score
We create a graph using the free online databases DisGeNET and DGIdb. The user must execute the R code with one argument (disease name as found in the DisGeNET_xxxxxxxx.txt) and a score-sorted list of potential related drugs through the common disease-drug genes is returned.

#Requires:
library(igraph)
DisGeNET_xxxxxxxx.txt and DGIdb_xxxxxxxx.txt in the same folder with the R script.

#Example:
Rscript disease-drug-score.R "Idiopathic Pulmonary Fibrosis"

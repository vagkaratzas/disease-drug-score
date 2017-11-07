#!/usr/bin/env Rscript

### Libraries
suppressMessages(library(igraph))

### Main
### User Inputs
args = commandArgs(trailingOnly = TRUE) #allow use of args
if (length(args) == 0) { #check if correct number of args
	stop("Need disease node string as argument, as written in DisGeNET.txt", call.=FALSE) #if not disease name, throw error and stop
}
disease <- args[1] #user input disease

DisGeNET <- read.delim("DisGeNET-07112017.txt", header = TRUE, sep = "\t") #read the DisGeNET file
DGIdb <- read.delim("DGIdb-07112017.txt", header = TRUE, sep = "\t") #read the DGIdb file
DisGeNET_edgelist <- DisGeNET[,c(1,2)] #keep gene-disease columns
DGIdb_edgelist <- DGIdb[,c(1,2)] #keep gene-disease columns
colnames(DisGeNET_edgelist) <- c(1,2) #need same colnames to rbind
colnames(DGIdb_edgelist) <- c(1,2) #need same colnames to rbind
edgeList <- rbind(DisGeNET_edgelist, DGIdb_edgelist) #combine to main edgelist
netNodes <- cbind(as.character(edgeList[,1]),as.character(edgeList[,2])) #make sure all nodes will be interpreted as strings
netNodes[is.na(netNodes)] <- "NA" #if value NA found, replace it with the string NA
graph <- graph_from_edgelist(netNodes, directed = FALSE) #make graph from massive edgelist
graph <- simplify(graph) #removes possible reoccuring edges, case sensitive

nodeIndex <- match(disease,V(graph)$name) #find index of given disease
neighbors <- as.matrix(neighbors(graph, nodeIndex)) #find the first level neighbor-genes of the disease

drug_list <- c() #initializing drug list
for (i in neighbors) { #for all neighbors of input disease check drug-neighbors
	gene_neighbors  <- as.matrix(neighbors(graph, i)) #find the first level neighbors of the current gene (will only keep drugs)
	#print(V(graph)[i]$name) #printing current gene name
	for (j in gene_neighbors) { #for all neighbors of current gene
		if (V(graph)[j]$name %in% DGIdb_edgelist$"2") { #only if neighbor is Drug
			#print(V(graph)[j]$name) #print name
			drug_list <- c(drug_list,V(graph)[j]$name) #add drug name to list
		}
	}
}

drug_times <- table(drug_list) #count genes that each drug targets
drug_scores <- drug_times/length(drug_times) #calculate score, hits/total genes targeted
print(drug_scores)

expfile <- paste(args[1], "drugScores.txt", sep = "_")
write.table(drug_scores, expfile, sep="\t")


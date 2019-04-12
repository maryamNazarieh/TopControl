###########################################
#  - TopControl  Project                  #
#  - September 2017, edited April 2018   #                   
#  - Copyright: Maryam Nazarieh           #
###########################################
setwd("~/Desktop/TopControl-master/Example/")
##########################################################
### input files #######
##########################################################
foldChange.file = "DESeq_res.csv"
networkw.File = "res_complete_TFmiR.csv"
mds.file = "TFmiRnw_mds.txt"
mcds.file = "TFmiRnw_mcds.txt"
hub.file = "TFmiRnw_hub.txt"

##########################################################
### output files ######
##########################################################
candidate.fourth.layer.file = "candidatesLayer4.csv"
candidate.fifth.layer.file = "candidatesLayer5.csv"
candidate.disease.file = "experimentallyValidatedNodes.csv"

##########################################################
### disease files ######
##########################################################
gene.disease.file = "gene_disease_DisGeNET.txt"
mirna.disease.file = "mir_disease_HMDD.txt"
disease.name = 'Breast Neoplasms'

##########################################################
### Node Degree Computation ######
##########################################################
library(igraph)
nw = read.delim(networkw.File,sep="\t")
graph = graph.data.frame(nw,directed=FALSE,vertices=NULL)
v <- get.data.frame(graph, what="vertices")
e <- get.data.frame(graph, what="edges")
d <- centralization.degree(graph,mode ="all")
d <- data.frame(gene=V(graph)$name, degree = d)
d <- d[order(-d$degree.res),]
##########################################################
### Read the input & disease files ######
##########################################################
nw_nodes <- v
#nw_nodes = read.table(network.nodes,sep="\t",col.names = FALSE) 
mds_nw = read.table(mds.file,sep="\t",col.names = FALSE)
mcds_nw = read.table(mcds.file,sep="\t",col.names = FALSE)
hub_nw = read.table(hub.file,sep="\t",col.names = FALSE)
DE_nodes = read.table(foldChange.file,sep="\t",header = TRUE)

gene.disease.table = read.table(gene.disease.file,sep="\t",quote="")
disease.gene.specific = gene.disease.table[grep(disease.name,gene.disease.table$V2,ignore.case=TRUE),]
mirna.disease.table = read.table(mirna.disease.file,sep="\t",quote="")
mirna.disease.specific = mirna.disease.table[grep(disease.name,mirna.disease.table$V3,ignore.case=TRUE),]
gene.mirna.disease.set = as.data.frame(c(as.character(disease.gene.specific[,1]),as.character(mirna.disease.specific[,2])))
gene.mirna.diseaseSet.nwNodes = intersect(nw_nodes[,1], gene.mirna.disease.set[,1])
##########################################################
### Read the input & disease files ######
##########################################################
score = rep(0,dim(nw_nodes)[1])
mds = rep(0,dim(nw_nodes)[1])
mcds = rep(0,dim(nw_nodes)[1])
hub = rep(0,dim(nw_nodes)[1])
disease_associated = rep(0,dim(nw_nodes)[1])
degree = rep(0,dim(nw_nodes)[1])
LFC = rep(0,dim(nw_nodes)[1])

##########################################################
### Score Calculation for Network Nodes using TopControl  
##########################################################
for (i in 1: dim(nw_nodes)[1]){
  if (nw_nodes[i,1] %in% mds_nw[,1]){
    score[i] = score[i] + 1
    mds[i] = mds[i] + 1
  }
  if (nw_nodes[i,1] %in% mcds_nw[,1]){
    score[i] = score[i] + 1
    mcds[i] = mcds[i] + 1
  }
  if (nw_nodes[i,1] %in% hub_nw[,1]){
    score[i] = score[i] + 1
    hub[i] = hub[i] + 1
  }
  if (nw_nodes[i,1] %in% gene.mirna.diseaseSet.nwNodes){
    disease_associated[i] = disease_associated[i] + 1
  }
  if (nw_nodes[i,1] %in% d$gene){
    res = d[which(d$gene %in% nw_nodes[i,1]),]$degree.res
    degree[i] = res
  }
  if (nw_nodes[i,1] %in% DE_nodes$id){
    LFC[i] = DE_nodes[which(DE_nodes$id %in% nw_nodes[i,1]),]$log2FoldChange
    LFC[i] = abs(round(LFC[i],digits = 2))
  }
  i = i + 1
}
gene_rank = cbind(nw_nodes,degree,hub,mds,mcds,score,LFC,disease_associated)
colnames(gene_rank) = c("gene","degree","hub_degree","mds","mcds","score","log2FoldChange","disease")

##########################################################
### TopControl Candidates in the fourth Layer 
##########################################################
candidates.layer.4 = cbind(nw_nodes,degree,hub,mds,mcds,score,LFC)
colnames(candidates.layer.4) = c("gene","degree","hub_degree","mds","mcds","score","log2FoldChange")
candidates.layer.4 <- candidates.layer.4[order(-score,-LFC),] 
candidates.layer.4 <- candidates.layer.4[which(candidates.layer.4$score >= 1),]
print(c("The number of Candidates in the fourth layer is:",dim(candidates.layer.4)[1]))
write.table(candidates.layer.4,candidate.fourth.layer.file,sep="\t",col.names = T,row.names = F)
##########################################################
### TopControl Candidates in the fifth Layer 
##########################################################
candidates.layer.5 <- gene_rank[order(-score,-LFC,-disease_associated),] 
candidates.layer.5 <- candidates.layer.5[which(candidates.layer.5$score == 3),]
print(c("The number of Candidates in the fifth layer is:",dim(candidates.layer.5)[1]))
write.table(candidates.layer.5,candidate.fifth.layer.file,sep="\t",col.names = T,row.names = F)

####################################################################
### TopControl-assigned scores of disease-associated genes and miRNAs
####################################################################
exp.validated.nodes <- gene_rank[order(-disease_associated,-score,-LFC),] 
name = noquote(disease.name)
exp.validated.nodes <- exp.validated.nodes[exp.validated.nodes$disease==1, ]
print(c("The number of experimentally validated Candidates is:",dim(exp.validated.nodes)[1]))
write.table(exp.validated.nodes,candidate.disease.file,sep="\t",col.names = T,row.names = F)

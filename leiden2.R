source('modelfns.R')
source('descend.R')
source('clustplots.R')
library(igraph)
library(dirfns)
library(fgsea)

groups <- read.csv('out/groups.csv')



int <- read.csv('interactions.csv',row.names=1)[,3:4]
g <- graph_from_edgelist(as.matrix(int),F)
g <- graph_from_adjacency_matrix(as_adjacency_matrix(g),'plus',T,diag=F)
graph.pdf('interactions',g)

write.graph(g,format='dot',file='tmp.dot')
out <- mkdate('interactions','dot')
system2('sed',paste("'s/name/label/' tmp.dot >",out))
system2('dot',paste('-Tsvg -O',out))

int <- rbind(int,setNames(int[,2:1],names(int)))
int <- rbind(as.matrix(int),t(sapply(unique(groups$Condition),rep,2)))
int <- int[!duplicated(int),]
dir.csv(int,'interactions.csv')

int <- apply(int,1,paste,collapse='->')

clusts2 <- read.clusts('2022-10-03/encode2')

k <- get.k(3:53,clusts2$dists,groups$Condition,int)
dir.csv(k,'k')
k.optim <- k[which.max(k[,2]),1]

knn <- get.knn(clusts$dists,k.optim,'plus')
dir.pdf('interaction_enrichment')
plotEnrichment(int,table(edge.factor(group.edge(knn,groups$Condition),F)))
dev.off()

leidens <- res.unif(c(0.01,3),k.optim,clusts2$encoded,int,groups$Condition,1000)
dir.csv(leidens,'leiden')

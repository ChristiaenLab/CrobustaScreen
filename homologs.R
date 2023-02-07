homologs <- read.csv("mmusculusHomologs.csv")
interactions <- read.csv('interactions.csv')
row.names(homologs) <- homologs$protein_external_id

interactions$from <- homologs[interactions$from,'name']
interactions$to <- homologs[interactions$to,'name']
interactions <- split(interactions,apply(interactions[,1:2],1,paste,collapse='_'))
interactions <- sapply(interactions,function(x) x[which.max(x[,3]),3])

interactions <- cbind(as.data.frame(do.call(rbind,strsplit(names(interactions),'_'))),interactions)

mat <- sapply(colnames(pois$log2OR),function(x) sapply(row.names(pois$log2OR),function(y) do.call(function(z){if(length(z)==0) 0 else z},list(interactions[interactions[,1]==x&interactions[,2]==y,3]))))

lfc <- log2(mat/mean(mat))
lfc[!is.finite(lfc)] <- 0

g <- graph_from_adjacency_matrix(lfc,'directed',T,F)

networkPois(g,'homologInteractions',title='log2FC',colfn = col.abs(lfc,0.01,c('white','red')))


odds <- pois$log2OR
odds[!is.finite(odds)] <- 0

dir.pdf('interactionRatio',opts$clusts,append.date = F)
Heatmap(odds-lfc,cluster_rows = F,cluster_columns = F)
dev.off()

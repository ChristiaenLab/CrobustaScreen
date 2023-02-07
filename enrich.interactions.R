gsea.prot <- function(k,dists,groups,int){
	require(fgsea)
	require(igraph)

	g <- graph_from_adjacency_matrix(get.knn(dists,k),'directed',T,F)

	ge <- get.edgelist(g)
	edges <- cbind(groups[ge[,1]],groups[ge[,2]])
	edges <- apply(edges,1,paste,collapse='->')
	ct <- fgsea(list(edge.ct=int),table(edges))
	return(ct$ES)
}


gsea.edges <- function(k,clusts,dists,groups){
	require(fgsea)
	require(igraph)

	g <- graph_from_adjacency_matrix(get.knn(dists,k),'directed',T,F)
	int <- read.csv('interactions.csv',row.names=1)
	int <- apply(int[,3:4],1,paste,collapse='->')

	ge <- get.edgelist(g)
	edges <- cbind(groups[ge[,1]],groups[ge[,2]])
	edges <- apply(edges,1,paste,collapse='->')
	gw <- get.edge.attribute(g)$weight
	names(gw) <- edges

	res <- fgsea(list(edges=int),gw)
	ct <- fgsea(list(edge.ct=int),table(edges))

	cl <- mapply(gsea.clust, colnames(clusts),
		     MoreArgs = list(clusts,edges,ge,int),SIMPLIFY = F)

	return(rbind(res,ct,do.call(rbind,cl)))
}

gsea.clust <- function(name,clusts,edges,nodes,int){
	edges <- as.factor(edges)
	clust <- clusts[,name]
	clust.from <- clust[nodes[,1]]
	clust.to <- clust[nodes[,2]]

	is.cis <- mapply('==',as.data.frame(clust.from),as.data.frame(clust.to))
	scores <- table(edges[is.cis])/table(edges[!is.cis])
	scores[!is.finite(scores)] <- max(Filter(is.finite,scores))
	cl <- fgsea(setNames(list(int),name),scores)
}

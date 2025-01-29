get.knn <- function(dists,k,mode='directed'){
	require(igraph)
	neighbors <- apply(dists,2,order)
	adj <- sapply(1:ncol(neighbors),function(i){
			      r <- dists[,i]
			      sel <- neighbors[-1:-k,i]
			      r[sel] <- 0
			      return(r)
	})
	g <- graph_from_adjacency_matrix(adj,mode,T,F)
}

edge.factor <- function(e,as.factor=T) {
	res <- apply(e,1,paste, collapse='->')
	if(as.factor) res <- as.factor(res)
	return(res)
}

group.edge <- function(g,groups){
	require(igraph)
	e <- as_edgelist(g)
	e <- cbind(groups[e[,1]],groups[e[,2]])
	return(e)
}

graph.pdf <- function(out, g,
		      layout = layout_nicely,
		      ..., append.date=T){
	dir.pdf(out, append.date = append.date)
	plot(g, layout = layout, ...)
	dev.off()
}



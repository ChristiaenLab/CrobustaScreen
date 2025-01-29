source("R/test.knn.R")

test.leiden <- function(res,k,g,dat,reps,groups,int,...){
	require(cluster)
	require(leiden)

	dists <- as.matrix(dist(dat))

	clust <- leiden(g,resolution_parameter=res)
	if(all(clust==1)) return(c(rep(0,6),clust))

	err <- test.knn(dat,k,clust,reps=reps)

	es <- test.leiden.gsea(g,clust,groups,int)
	sil <- mean(silhouette(clust,dists)[,3])

	g.gene <- gene.network(g,res,groups)
	recall <- score.gene.network(g.gene,int)

	#stat <- err*es*sil*recall
	stat <- es * recall
	return(c(combined_score = stat,
		 ES = es, 
		 recall = recall,
		 log2error = err,
		 mean_silhouette = sil,
		 nclust = max(clust), 
		 clust))
}

test.leiden.gsea <- function(g,clust,groups,int,...){
	require(fgsea)
	require(leiden)

	e.group <- edge.factor(group.edge(g,groups),T)

	e.clust <- group.edge(g,clust)

	is.cis <- e.clust[,1]==e.clust[,2]
	scores <- table(e.group[is.cis])/table(e.group[!is.cis])
	scores[!is.finite(scores)] <- max(Filter(is.finite,
						 scores))
	cl <- fgsea(list(interactions=int),scores)

	return(cl$ES)
}

gene.network <- function(g, res, groups,
			 mode = 'directed', 
			 weighted = T,
			 diag = F, ...){
	require(igraph)
	gene.edge <- group.edge(g, groups)
	K.out <- sapply(split(igraph::degree(g, mode='out'),
			      groups), sum)
	K.in <- sapply(split(igraph::degree(g, mode='in'),
			     groups), sum)

	fy <- function(y, x){ 
		e <- sum(gene.edge[,1] == x & gene.edge[,2] == y) 
		e - (res * ((K.out[[x]] * K.in[[y]]) / 
			    (2 * sum(igraph::degree(g)))))
	}

	fx <- function(x) { 
		sapply(unique(groups), fy, x)
	}

	clust.modularity <- sapply(unique(groups), fx)

	clust.modularity[clust.modularity < 0] <- 0
	g.gene <- graph_from_adjacency_matrix(clust.modularity,
					      mode, weighted,
					      diag, ...)
	return(g.gene)
}



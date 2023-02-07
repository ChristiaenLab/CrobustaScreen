get.knn <- function(dat,k){
	require(igraph)
	dists <- as.matrix(dist(dat))
	neighbors <- apply(dists,2,order)
	adj <- sapply(1:ncol(neighbors),function(i){
			      r <- dists[,i]
			      sel <- neighbors[-1:-k,i]
			      r[sel] <- 0
			      return(r)
	})
	g <- graph_from_adjacency_matrix(adj,'directed',T,F)
	return(g)
}

get.clusts <- function(dat,k,res=0.5,...){
	require(leiden)
	g <- get.knn(dat,k)
	clust <- leiden(g,resolution_parameter=res,...)
	return(clust)
}


test.knn <- function(dat,k,clust,reps=1000,...){
	trainsel <- sample(1:nrow(dat),nrow(dat) %/% 10 * 9)
	traindat <- dat[trainsel,]
	testdat <- dat[-trainsel,]
	cv <- replicate(reps,knn.err(dat,k,clust))
	cv[!is.finite(cv)] <- -log2(1/(2*nrow(testdat)))
	err <- mean(cv)
	return(err)
}

knn.err <- function(dat,k,clust){
	require(class)

	trainsel <- sample(1:nrow(dat),nrow(dat) %/% 10 * 9)
	traindat <- dat[trainsel,]
	testdat <- dat[-trainsel,]

	cv <- class::knn(traindat,testdat,clust[trainsel],k)
	err <- sum(cv!=clust[-trainsel])/nrow(testdat)
	return(-log2(err))
}

source("R/descend.R")
source("R/test.knn.R")
source("R/leiden.R")

par.apply <- function(..., f = pbsapply){
	require(parallel)
	require(pbapply)

	ncore <- detectCores() - 2
	cl <- makeCluster(ncore, "FORK")
	res <- f(cl = cl, ...)
	stopCluster(cl)
	return(res)
}

# I got way too clever with this
eta <- function(f, pb){
	require(progress)
	function(...){
		pb$tick()
		f(...)
	}
}

map.eta <- function(FMAP, msg){ 
	require(progress)
	require(stringr)

	fmt = str_interp("  ${msg} [:bar] :percent eta: :eta")
	function(FUN, X, ...){
		total <- length(X)
		pb <- progress_bar$new(format = fmt, 
				       total = total, 
				       clear = FALSE, 
				       width= 60)
		f  <- eta(FUN, pb)
		FMAP(FUN = f, X, ...)
	}
}

get.res <- function(range,breaks,rate,maxiter,
		    k,dat,int,groups,...){
	require(igraph)
	dists <- as.matrix(dist(dat))
	g <- get.knn(dists,k)
	descend(range,test.leiden,
		k=k,g=g,dat=dat,reps=1000,groups=groups,
		int=int,breaks=breaks,rate=rate,
		maxiter=maxiter,min=0.05,
		trace=matrix(nrow=0,ncol=7+nrow(dat)),...)
}

get.k <- function(range, dists, group, int, mode='plus'){
	res <- par.apply(range, gsea.prot, 
			 dists, group, int, 
			 mode = mode)
	return(cbind(k = range, ES = res))
}

get.res.unif <- function(range, k, dat, int, groups,
			 reps = 1000){
	require(parallel)
	require(igraph)

	dists <- as.matrix(dist(dat))
	g <- get.knn(dists, k, 'plus')

	res <- runif(reps, range[1], range[2])
	out <- par.apply(res, test.leiden,
		k = k,g = g,dat = dat,reps = 1000,
		groups = groups,
		int = int)
	return(cbind(resolution = res, t(out)))
}

get.leidens <- function(range, k, dat){
	require(igraph)
	require(leiden)
	 
	dists <- as.matrix(dist(dat))
	g <- get.knn(dists,k,'plus')

	res <- runif(reps,range[1],range[2])
	clust <- par.apply(leiden,resolution_parameter=res,
			   MoreArgs=list(object=g),
			   f=clusterMap)

	return(cbind(res,t(clusts)))
 }


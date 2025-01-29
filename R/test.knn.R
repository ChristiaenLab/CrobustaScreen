source("R/hyper.R")
source("R/knn.R")

test.knn <- function(dat,k,clust,sample=0.9,reps=1000,...){
	n.train <- round(nrow(dat)*sample)
	cv <- replicate(reps,knn.err(dat,k,clust,n.train))
	cv[!is.finite(cv)] <- 1/(2+n.train)
	err <- mean(cv)
	mse <- mean(cv^2)
	#         return(c(error=err,MSE=mse))
	return(-log2(err))
}

knn.err <- function(dat,k,clust,n.train){
	require(class)

	trainsel <- sample(1:nrow(dat),n.train)
	traindat <- dat[trainsel,]
	testdat <- dat[-trainsel,]

	cv <- class::knn(traindat,testdat,clust[trainsel],k)
	err <- sum(cv!=clust[-trainsel])/n.train
	return(err)
}

graph.hyper <- function(g,conds,x,y){
	require(igraph)
	deg <- igraph::degree(g)
	gene.edge <- group.edge(g,conds)
	e <- sum(gene.edge[,1]==x&gene.edge[,2]==y) 
	m <- sum(deg[conds==x])
	k <- sum(deg[conds==y])
	n <- sum(deg)-m
	p <- phyper(e-1,m,n,k)
	odds <- (e/(k-e))/(m/n)
	return(c(log2OR=odds,p=p,q=e))
}

get.hyper <- function(g,conds,padj.method="fdr"){
	require(purrr)
	condl <- unique(conds)
	hyper <- mapply(partial(mapply,
				partial(graph.hyper,knn,
					conds),
				condl),
			condl,SIMPLIFY=F)
	hyper <- sapply(c('log2OR','p','q'),
			partial(sapply,hyper,'['),,
			simplify=F)
	hyper$p <- apply(hyper$p,2,p.adjust,
			 method=padj.method,
			 n=length(condl))
	return(hyper)
}

score.gene.network <- function(g,int){
	e.gene <- edge.factor(as_edgelist(g),F)
	recall <- length(intersect(e.gene,int))/length(int)
	return(recall)
}

gsea.prot <- function(k,dists,groups,int,mode='directed'){
	require(fgsea)
	require(igraph)

	g <- get.knn(dists,k,mode=mode)

	ge <- as_edgelist(g)
	edges <- cbind(groups[ge[,1]],groups[ge[,2]])
	edges <- apply(edges,1,paste,collapse='->')
	ct <- fgsea(list(edge.ct=int),table(edges))
	return(ct$ES)
}


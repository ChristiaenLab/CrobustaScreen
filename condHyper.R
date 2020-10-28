#' Hypergeometric test for enrichment of conditions in a cluster.
#' @param id A vector of sample IDs.
#' @param conds A vector of the same length as ID giving the condition of each sample.
#' @param clusts A vector of the same length as ID giving the cluster ID of each sample.
#' @param padj.method Method passed to \code{\link{p.adjust}} for multiple hypothesis correction.
#' @seealso \code{\link{phyper2}}, \code{\link{p.adjust}}
#' @export
condHyper <- function(id,conds,clusts,padj.method='fdr'){
	#         id <- id[!is.na(conds)]
	#         cols <- unique(clusts)
	#         clusts <- clusts[!is.na(conds)]
	test <- split(id,conds)
	clusts <- split(id,clusts)
	#         fn <- function(x) sum(!is.na(x))
	m <- sapply(test,length)
	n <- length(id)-m
	k <- sapply(clusts,length)
	q <- as.data.frame(sapply(clusts,function(m) sapply(test, function(k){
		sum(m%in%k)
	})))
	log2OR <- mapply(
	  function(q.i,k.i) mapply(
	    function(q.ij,m.j) log2(
	      (q.ij/(k.i-q.ij))/(m.j/(length(id)-m.j))
	    ),
	    q.i,m
	  ),
	  q,k
	)
	row.names(log2OR) <- names(test)

	testHyper <- mapply(function(q,k) mapply(
	  phyper2,q=q-1,k=k,m=m,n=n
	),q=q,k=k)
	testFdr <- apply(testHyper,2,p.adjust,method=padj.method)
	row.names(testHyper) <- names(test)
	return(list(log2OR=log2OR,FDR=testFdr,q=q))
}

#' Two-tailed version of \code{\link{phyper}}. It simply takes the minimum of the upper and lower tails and multiplies the result by 2.
#' @param ... Arguments to \code{phyper()}.
#' @export
#' @seealso \code{\link{phyper}}
phyper2 <- function(...) min(
	phyper(...,lower.tail=T),
	phyper(...,lower.tail=F)
)*2


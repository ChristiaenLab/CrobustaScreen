getPois <- function(cond, dists){
	#select distances between points in the same condition
	condsel <- sapply(unique(cond), function(x) cond == x)
	#condition permutations to test
	conds <- combn(unique(cond), 2)
	conds <- cbind(sapply(unique(cond), rep, 2), conds)

	#number of embryos per condition
	ncond <- table(cond)
	#neighbors have dist > 0
	n <- dists > 0

	#get matrix of neighbors for each cond
	kn <- split(as.data.frame(dists > 0), cond)

	#poisson test
	pois <- lapply(kn, function(x){
		       #split columns by condition
		       conds <- split(as.data.frame(t(x)), cond)
		       #total number of neighbors
		       k <- sum(as.matrix(x))
		       #test each condition
		       res <- lapply(conds, function(y){
				y <- as.matrix(y)
				#number of neighbors in test condition
				ct <- sum(y)
				#proportion of test condition to total embryos
				p <- nrow(y)/sum(ncond)
				pval <- poisson.test(ct, k, p)$p.value
				#proportional overrepresentation of test condition among all neighbors vs. expected
				odds <- log2((ct/k)/p)
				return(c(p = pval, log2OR = odds, count = ct))
			})
		       return(do.call(rbind, res))

		})
	odds <- sapply(pois, '[', ,'log2OR')
	fdr <- sapply(pois, '[', ,'p')
	qval <- sapply(pois, '[', ,'count')
	return(list(log2OR = odds, p = fdr, count = qval))
}



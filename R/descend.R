# Approximates gradient descent for nonanalytic functions.
# Samples from a given `range`. After each iteration,
# `range` is recentered on the optimal value
# and shrunk by learning rate `rate`.

descend <- function(range,f,...,breaks,rate,maxiter,
		    x=0,fwhich=which.max,discrete=F,min=0,
		    trace=matrix(nrow=0,ncol=2),parallel=T){
	require(parallel)
	if(parallel){
		applyfn <- function(...){
			ncore <- min(breaks,detectCores()-2)
			cl <- makeCluster(ncore,"FORK")
			res <- parSapply(cl,...)
			stopCluster(cl)
			return(res)
		}
	}else{applyfn <- sapply}

	res <- list(argmax=x,trace=trace)
	if(nrow(trace)/breaks==maxiter) return(res)

	if(discrete){
		i <- setdiff(round(range[1]):round(range[2]),
			     trace[,1])
		if(length(i)==0) return(res)
		if(length(i)>breaks){
			i <- sample(i,breaks)
		}
	}else{
		i <- runif(breaks,range[1],range[2])
	}
	print(as.character(i))

	vals <- applyfn(i,f,...)
	if(length(dim(vals)>1)) vals <- t(vals)
	res$trace <- rbind(trace,cbind(i,vals))
	res$trace <- res$trace[order(res$trace[,1]),]
	sel <- fwhich(res$trace[,2])
	res$argmax <- res$trace[sel,1]
	win <- (range[2]-range[1])/2*rate

	range.o <- c(res$argmax-win,res$argmax+win)
	if(range.o[1]<min) range.o[1] <- min

	descend(range.o,f,...,breaks=breaks,rate=rate,
		maxiter=maxiter,x=res$argmax,fwhich=fwhich,
		discrete=discrete,min=min,trace=res$trace)
}


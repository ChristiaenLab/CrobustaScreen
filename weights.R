library(keras)
library(functional)
library(parallel)
library(ggplot2)
library(circlize)

dat <- read.csv('out/z_dat.csv',row.names=1)
dat <- dat[,apply(abs(dat),2,max)!=0]

encode2 <- load_model_hdf5('2023-01-23/encode2.model')

outmodel <- keras_model(inputs=get_layer(encode2,
					 index=7)$input,
			outputs=get_layer(encode2,
					  index=12)$output)

predict(outmodel,diag(1,2,2))
embedding <- do.call(rbind,lapply(runif(100,-1,1),cbind,runif(100,-1,1)))
pred <- predict(outmodel,embedding)

dat <- setNames(as.data.frame(cbind(embedding,pred)),
		  c('embedding1','embedding2',names(dat)))
dat <- dat[seq(10,10000,10),]

plots <- lapply(names(dat)[-1:-2],
		function(y) ggplot(dat,
				   aes_string(x=names(dat)[1],
					      y=y,
					      col=names(dat)[2]))+
		scale_color_viridis_c()+geom_point())

arrange.stats(plots,'embeddings',height=128,width=16)
arrange.stats(plots[1:12],'embeddings',height=12,width=16)

library(rgl)

get.nucpos <- function(x) 
	do.call(cbind, 
		split(t(do.call(rbind,
				lapply(c('+','-'),
				       mapply, 
				       x[,81:83], 
				       x[,66:74]))) ,1:3))
get.cellpos <- function(x)
	do.call(cbind, 
		split(t( x[,24:32]) ,1:3))


plot3d(ellipse3d(get.nucpos(celldat[[1]])))
surface3d(get.nucpos(celldat[[1]][1:5,]),)
rglwidget()

weight <- encode2$get_weights()
weight <- weight[sapply(weight,Compose(dim,length))==2]

fargs <- function(f,...,argn='') {
	argsx <- list(...)
	function(...) {
		argsy <- list(...)
		if(argn!='') names(argsy)[1] <- argn
		args <- append(argsx,argsy)
		do.call(f,args)
	}
}

outmodel <- do.call(keras_model,lapply(6:12,function(i) get_layer(encode2,index=i)))#fargs(get_layer,encode2,argn='index')))

lsiter <- function(l,f,...,cond=is.list){
	if(cond(l)) lapply(l,lsiter,f,...)
	else f(l,...)
}


ncore <- detectCores()-2
cl <- makeCluster(ncore,"FORK")

get.paths <- function(cl,x,y){
	if(length(dim(x))==2){
	parLapply(cl,1:ncol(x),function(i) t(sapply(x[,i],'*',y[i,])))
	}
	else lapply(x,get.paths,y,cl=cl)
}

paths <- Reduce(fargs(get.paths,cl),weight[3:5],weight[[2]])

stopCluster(cl)
fargs(get.paths,cl)(weight[[1]],weight[[2]])

lsiter(paths,function(x) x[1,])

contrib <- lapply(1:114,function(i) Reduce(function(x,y){
		if(length(dim(x)

paths <- Reduce(function(x,y){
vals <- matrix(1,114,114)
diag(vals) <- 0

predict(encode2,diag(nrow=114,ncol=114))

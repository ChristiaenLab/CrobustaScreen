source('descend.R')
source('plotfns.R')

library(leiden)
library(optparse)
library(igraph)
library(dirfns)
library(purrr)

parser <- OptionParser()
parser <- add_option(parser, '--modeldir', action = 'store',
		     default = '2023-02-01')
parser <- add_option(parser, '--clustdir', action = 'store',
		     default = '2023-02-01')
parser <- add_option(parser, '--groups', action = 'store',
		     default = 'out/groups.csv')
parser <- add_option(parser, '--pheno', action = 'store',
		     default = 'out/phenotype.csv')
opts <- parse_args(parser)


groups <- read.csv(opts$groups)
pheno <- read.csv(opts$pheno,row.names=1)

dat <- as.matrix(read.csv(paste0(opts$modeldir,"/dat.csv"),row.names=1))

model <- read.model(paste0(opts$modeldir,'/encode2.model'),dat)

model$k <- read.csv(paste0(opts$clustdir,'/k.csv'),row.names=1)
names(model$k)[1:2] <- c('k','ES')


clusts <- lrtab(opts$clustdir,read.csv,'leiden',row.names=1)

polycv <- function(x,y,d){
	train(poly(x,d),y,method="lm",
	      trControl=trainControl(method="LOOCV"))
}

res.fit <- lapply(1:20,partial(polycv,
			       clusts$leiden$res,
			       clusts$leiden$ES))

lapply(clusts,function(x){
	       train(ES ~ poly(res,1),
		     method="lm",data=x,
		     trControl=trainControl(method="LOOCV"))
		     })

# which.k <- which.max(sapply(clusts, 
#                       purrr::compose(max, '['), ,
#                       'stat'))
# k <- as.numeric(sub('leiden.k','',names(clusts)[which.k]))

which.k <- which.max(model$k$ES)
k <- model$k$k[which.k]

knn <- get.knn(model$dists,k,'plus')

if(length(clusts)>1){
	model$clusts <- clusts[[which.k]]
}else model$clusts <- clusts[[1]]
names(model$clusts)[1:2] <- c('resolution','combined_score')

model$clusts <- model$clusts[model$clusts$nclust>1,]

mapply(statplot,clusts,sub('leiden.','',names(clusts)))

plot.edge(model$encoded,knn,'knn')

dir.f(ggexport)(dot.col('encoding2',
			as.data.frame(model$encoded),
			col=groups$Condition,'condition'),
		filename='encoding.pdf')

plots <- lapply(names(model$clusts)[2:7],dot.stat,clusts$clusts)

es <- dot.stat('ES',model$k)
dir.f(ggexport)(es,filename="ES.pdf")


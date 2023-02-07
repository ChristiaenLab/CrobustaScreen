source('modelfns.R')
source('descend.R')
source('clustplots.R')

groups <- read.csv('out/groups.csv')
clusts <- read.model('2023-10-03/encode2')

clusts$clusts <- read.csv('2022-10-22/leiden.csv',row.names=1)
names(clusts$clusts)[1:2] <- c('resolution','combined_score')
clusts$clusts <- clusts$clusts[clusts$clusts$nclust>1,]

clusts$k <- read.csv('2022-10-22/k.csv',row.names=1)
names(clusts$k)[1:2] <- c('k','ES')

k <- clusts$k$k[which.max(clusts$k$ES)]
knn <- get.knn(clusts$dists,k,'plus')

err <- par.apply(test.knn,clust=as.data.frame(t(clusts$clusts[,-1:-7])),
		 MoreArgs=list(dat=clusts$encoded,k=k),f=clusterMap)

err <- do.call(rbind,err)

mod.pheno <- par.apply(modularity,
		       resolution=clusts$clusts[,1],
		       MoreArgs=list(x=knn,
				     membership=as.factor(groups$Phenotype)),
		       f=clusterMap)

mod <- par.apply(modularity, 
		 resolution=clusts$clusts[,1], 
		 membership=as.data.frame(t(clusts$clusts[,-1:-7])), 
		 MoreArgs=list(x=knn),
		 f=clusterMap)

dir.csv(cbind(clusts$clusts[,1:7],err,modularity=unlist(mod),phenotype.modularity=unlist(mod.pheno)),'stats')

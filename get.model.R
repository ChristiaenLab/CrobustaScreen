source("imageFns.R")
source('modelfns.R')
library(dirfns)
library(moreComplexHeatmap)

groups <- read.csv('out/groups.csv')
dat <- read.csv('out/z_dat.csv',row.names=1)
dat <- cbind(dat,zdf(groups[,c("Threshold.nuclei","Threshold.membrane")]))

dat <- dat[,apply(abs(dat),2,max)!=0]
x_train <- t(t(dat)/apply(abs(dat),2,max))
cols <- groups[,c("Condition","Phenotype")]

nuc.sel <- which(x_train[,115]==0)
mem.sel <- which(x_train[,116]==0)
sel <- union(nuc.sel,mem.sel)

denoise.nuc <- denoise(x_train[-sel,-115],
		       x_train[-sel,115],
		       c(115,58,29,14,8,4,2,1),
		       'denoise.nuc')
denoise.mem <- denoise(x_train[-sel,-116],
		       x_train[-sel,116],
		       c(115,58,29,14,8,4,2,1),
		       'denoise.mem')
denoise.thr <- denoise(x_train[-sel,-115:-116], 
		       x_train[-sel,115:116], 
		       c(114,58,29,14,8,4,2), 
		       'denoise.thr')

x_train[mem.sel[1],116] <- predict(denoise.mem,x_train[mem.sel[1],-116,drop=F])
x_train[nuc.sel,115:116] <- predict(denoise.thr,x_train[nuc.sel,-115:-116,drop=F])

dir.csv(x_train,'dat')
		
encode2 <- autoencode(x_train,c(58,29,14,8,4,2),
		      cols,'encode2',epochs=10000)
encode3 <- autoencode(x_train,c(58,29,14,7,3),
		      cols,'encode3',epochs=10000)
encode7 <- autoencode(x_train,c(58,29,14,7),
		      cols,'encode7',epochs=10000)
encode14 <- autoencode(x_train,c(58,29,14),
		      cols,'encode14',epochs=10000)


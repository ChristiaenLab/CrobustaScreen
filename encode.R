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

# nuc.sel <- which(x_train[,115]==0)
# mem.sel <- which(x_train[,116]==0)
# sel <- union(nuc.sel,mem.sel)
# 
# denoise.nuc <- denoise(x_train[-sel,-115],
#                        x_train[-sel,115],
#                        c(115,58,29,14,8,4,2,1),
#                        'denoise.nuc')
# denoise.mem <- denoise(x_train[-sel,-116],
#                        x_train[-sel,116],
#                        c(115,58,29,14,8,4,2,1),
#                        'denoise.mem')
# denoise.thr <- denoise(x_train[-sel,-115:-116], 
#                        x_train[-sel,115:116], 
#                        c(114,58,29,14,8,4,2), 
#                        'denoise.thr')
# 
# x_train[mem.sel[1],116] <- predict(denoise.mem,x_train[mem.sel[1],-116,drop=F])
# x_train[nuc.sel,115:116] <- predict(denoise.thr,x_train[nuc.sel,-115:-116,drop=F])

dir.csv(x_train,'dat')

hm <- Heatmap(t(x_train),name='standardized z-score',
	     #              cell.h=0.024, cell.w=0.005, 
	     column_split=cols$Condition,
	     show_column_names=F,
	     height = unit(.8,"npc"),
	     width=unit(.7,"npc"),
	     column_title_rot = 90)

dir.pdf('params',height=20,width=20)
draw(hm)
dev.off()

encode2 <- autoencode(x_train,c(58,29,14,8,4,2),
		      cols,'encode2',epochs=10000)
encode3 <- autoencode(x_train,c(58,29,14,7,3),
		      cols,'encode3',epochs=10000)
encode7 <- autoencode(x_train,c(58,29,14,7),
		      cols,'encode7',epochs=10000)
encode14 <- autoencode(x_train,c(58,29,14),
		      cols,'encode14',epochs=10000)

models <- list(encode2,encode3,encode7,encode14)

err <- sapply(models,evaluate,x_train,x_train)

layers <- sapply(models,function(x) length(x$layers))

encoded <- mapply(function(x,y) model.out(x,x_train,y/2),models,layers)

aic <- mapply(function(k,L) 2*k-2*log(1-L),sapply(encoded,ncol),err)

tab <- data.frame(bottleneck=sapply(encoded,ncol),layers=layers,MSE=err,AIC=aic)

dir.csv(tab,'aic')



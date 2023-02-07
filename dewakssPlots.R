library(optparse)
source('clustplots.R')

parser <- OptionParser()
parser <- add_option(parser, '--params', action = 'store',default = 'out/params.csv')
parser <- add_option(parser, '--clusts', action = 'store',default = '2022_03_10/raw/')
opts <- parse_args(parser)

# data from DEWAKSS
clusts <- paste0(opts$clusts, '/dat/obs.csv')
clusts <- read.csv(clusts,row.names=1)
# input data
dat <- read.csv(opts$params, row.names=1)
dat <- dat[,sapply(dat,var)!=0]
dat <- dat[,-105:-106]

umap <- read.csv(paste0(opts$clusts, '/dat/obsm.csv'))
row.names(umap) <- row.names(clusts)
pcs <- umap[,-ncol(umap)+1:-ncol(umap)]

distmat <- as.matrix(read.csv(paste0(opts$clusts, '/dists.csv'), header=F))
row.names(distmat) <- row.names(clusts)
colnames(distmat) <- row.names(clusts)

test <- clusts[,grep('leiden', names(clusts))]
out <- paste0(opts$clusts, '/', names(test))

pois <- getPois(clusts$Condition,distmat)

clustpcs <- lapply(test,function(x){
	sapply(split(pcs,x),sapply,mean)
})

library(moreComplexHeatmap)

quantHeatmap(pcs, 'pcs', cell.w=0.005, cell.h=0.005, 
	     conds=clusts[,c(1,3:5)],path=opts$clusts,
	     show_row_names=F)

quantHeatmap(pcs, 'pcClust', cell.w=0.005, cell.h=0.005, 
	     conds=clusts,path=opts$clusts,
	     show_row_names=F)

mapply(function(x,y){
        quantHeatmap(x,paste0(y,'/clustPCs'),quant=0.05,cell.w=0.1,cell.h=0.2)
#	dir.pdf('clustPCs',y,F)
#	draw(Heatmap(log2(x)))
#	dev.off()
},clustpcs,out)

hm <- clustpcs[grep('leiden_denoised',names(clustpcs))]
sp <- mapply(rep,names(hm),sapply(hm,ncol))
hm <- do.call(cbind,hm)
quantHeatmap(hm, 'pcDenoised', cell.w=0.500, cell.h=0.500, 
	     column_split=unlist(sp),path=opts$clusts,
	     show_row_names=F)

plotDist(umap,clusts$Condition,distmat,opts$clusts)

enrichCond(clusts$Condition,distmat,opts$clusts)

mapply(condDist, test, out=out, MoreArgs=list(clusts$Condition, as.matrix(dist(pcs))))

sil <- mapply(silhouettePlot, clusts=test, out=out, MoreArgs=list(dat=pcs, distmat=distmat),SIMPLIFY = F)

mapply(clustplots, clusts=test, out=out, MoreArgs=list(dat=umap, conds=clusts$Condition, distmat=distmat))

mapply(clusthyper, clusts=test, out=out, MoreArgs=list(dat=clusts))

mapply(clustparam, clusts=test, out=out, MoreArgs=list(dat=dat))

mse <- read.csv(paste0(opts$clusts,'/MSE.csv'),header=F)
msepcs <- sapply(split(as.data.frame(t(mse[-1:-2,-1:-2])),unlist(mse[1,-1:-2])),t,simplify=F)
msek <- sapply(split(as.data.frame(t(mse[-1:-2,-1:-2])),unlist(mse[2,-1:-2])),t,simplify=F)

dir.pdf('mse',opts$out,append.date=F,width=10,height=20)
Heatmap(t(sapply(msepcs,apply,2,mean)))
dev.off()

dir.pdf('paramPCcor',opts$out,append.date=F,width=10,height=20)
Heatmap(sapply(pcs,function(x) sapply(dat[,-1],cor, x,method='spearman')))
dev.off()

dir.pdf('silhouette', opts$out, append.date=F)
boxplot(sapply(sil,'[',,'sil_width'),las=2)
dev.off()

dir.pdf('silMeanClust', opts$out, append.date=F)
boxplot(lapply(lapply(sil,summary),'[[','clus.avg.widths'),las=2)
dev.off()

dir.pdf('silClust',opts$out,append.date=F)
boxplot(formula=sil_width ~ cluster,as.data.frame(sil))
dev.off()


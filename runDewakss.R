library(ComplexHeatmap)
library(circlize)
library(moreComplexHeatmap)
library(dirfns)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, '--params', action = 'store',default = 'out/params.csv')
parser <- add_option(parser, '--clusts', action = 'store',default = 'out/2020-10-14/raw/dat/obs.csv')
parser <- add_option(parser,'--out',action='store',default='out/2020-10-14/raw')
opts <- parse_args(parser)

# This script contains the dotplot function.
# It sources hmfns.R, which contains helper functions for plotting heatmaps,
# and dirfns.R, which contains helper functions for organizing outputs into directories.
# source("clustfns.R")
source('condHyper.R')

# data from DEWAKSS
clusts <- read.csv(opts$clusts,row.names=1)
# input data
dat <- read.csv(opts$params, row.names=1)
dat <- dat[,sapply(dat,var)!=0]

# unused features
# dat <- dat[row.names(clusts),c(-21:-26,-111,-116,-121:-131)]

# split data by clusters
clustdat <- split(dat,clusts$leiden)

# run hypergeometric tests for enrichment of conditions and phenotypes
hyper <- lapply(clusts[,c('Condition','Phenotype')], function(x) condHyper(row.names(clusts),x,clusts$leiden))

# condHM(hyper$odds,hyper$fdr,hyper$q,'condition',opts$out,row_split=rowsplit,row_title_rot=0)

#condition <- mapply(function(cond,id) {
#			    #             path <- paste0(path,'/',id)
#    cond <- condtest(row.names(umap$data),cond,lclust)
#    #             mapply(dir.csv,cond,names(cond),path)
#    #             condHM(cond$log2OR,cond$FDR,'condition',path)
#    return(cond)
#},conds,names(conds),SIMPLIFY=F)

# extract fields from test & reformat as matrices
odds <- do.call(rbind,lapply(hyper,'[[','log2OR'))
fdr <- do.call(rbind,lapply(hyper,'[[','FDR'))
qval <- do.call(rbind,lapply(hyper,'[[','q'))

# split phenotype & condition into separate panels
rowsplit <- unlist(
	mapply(
	       function(x,y) rep(y,nrow(x$log2OR)),
	       hyper,
	       names(hyper)
	)
)

# write matrices to csv
dir.csv(cbind(
	parameter=rowsplit,
	condition=row.names(odds),
	as.data.frame(odds)
),'log2OR', opts$out, append.date=F)
dir.csv(cbind(
	parameter=rowsplit,
	condition=row.names(odds),
	as.data.frame(fdr)
),'FDR', opts$out, append.date=F)
dir.csv(cbind(
	parameter=rowsplit,
	condition=row.names(odds),
	as.data.frame(qval)
),'size', opts$out, append.date=F)

# write heatmap
hyperDot(
	odds, 
	fdr, 
	qval, 
	file='condition', 
	path=opts$out, 
	row_split=rowsplit, 
	row_title_rot=0
)

dotPscale(
	odds, 
	fdr, 
	qval, 
	file='pscale', 
	path=opts$out, 
	row_split=rowsplit, 
	row_title_rot=0
)

# test for significance of each feature in each cluster
testfn <-function(x, clust, dat){
      mu <- mean(dat[,x],na.rm=T)
      FC <- mean(clust[,x],na.rm=T)/mu
      if(FC!=1&length(unique(clust[,x]))>1){
	      utest <- wilcox.test(clust[,x],mu=mean(dat[,x], na.rm=T))$p.value
	      ttest <- t.test(clust[,x],mu=mean(dat[,x], na.rm=T))$p.value
      }else{
	      utest <- NaN
	      ttest <- NaN
      }
      return(c(FC=FC,u=utest,t=ttest))
} 
clusttest <- function(clust){
       sapply(colnames(clust), testfn, clust, dat)
}
test <- lapply(clustdat, clusttest)

# select fields from output
fcdat <- sapply(test,'[',"FC",)
p.t <- sapply(test,'[','t',)
p.u <- sapply(test,'[','u',)
# fcdat <- sapply(test,sapply,'[[','FC')
# fcdat <- do.call(rbind,apply(fcdat,1,unlist))
# p.t <- sapply(test,sapply,'[[','t.p.value')
# p.u <- sapply(test,sapply,'[[','u.p.value')

# convert p-values to FDR values
fdr.t <- apply(p.t,2,function(x) p.adjust(unlist(x)))

# define color scale
clustdat <- lapply(clustdat, function(x) t(x)[row.names(fcdat),])
# col.z creates a color scale centered on 0 ranging from the 0.01 to 0.99 quantiles of the input data
# This scale will be used for the fold changes between average feature values in each cluster from the background.
colfc <- col.z(fcdat,.05,1)
# creates a color scale from 0 to 2
# This will give a log10 scale for FDR values between 1 and 0.01
colfdr <- colorRamp2(c(0,2),c('white','black'))

# show feature values within clusters as boxplots color-coded by FDR value
getAnn <- function(x) {
	m <- clustdat[[x]]
	fc.cols <- colfc(fcdat[,x])
	fdr.cols <- colfdr(fdr.t[,x])
	return(anno_boxplot(
		m, 
		which='row',
		#                 type='violin', 
		width=unit(1,'in'),
		box_width=0.9,
		#                 outline=F,
		gp=gpar(
			fill=fc.cols,
			col=fdr.cols
		)
	))
}

# apply annotation function to each cluster
ha <- lapply(colnames(fcdat),getAnn)
names(ha) <- paste('cluster',colnames(fcdat))
# bind annotation into single object
ha <- do.call(HeatmapAnnotation,append(ha,list(which='row')))

# assign names to color scales
lgd <- list(
	Legend(col_fun = colfc, title = "FC"),
	Legend(col_fun = colfdr, title = "-log10(FDR)")
)

hm <- hm.cell(fcdat,
	right_annotation = ha,
	cell.w=.15,
	cell.h=.15,
	show_heatmap_legend=F,
	col=colfc
)

dir.pdf('t', opts$out, height=24, width=24, append.date=F)
draw(hm, annotation_legend_list=lgd)
dev.off()


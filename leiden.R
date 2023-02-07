source('modelfns.R')
source('descend.R')
source('clustplots.R')

library(optparse)
library(igraph)
library(dirfns)

parser <- OptionParser()
parser <- add_option(parser, '--dir', action = 'store',
		     default = '2023-01-27')
parser <- add_option(parser, '--dat', action = 'store',
		     default = '2023-01-27/dat.csv')
parser <- add_option(parser, '--leiden_reps', 
		     action = 'store', default = 1000)
parser <- add_option(parser,'--all_k',action='store_true',
		     default=F)
opts <- parse_args(parser)

dat <- as.matrix(read.csv(opts$dat,row.names=1))

models <- lrtab(opts$dir,
		read.model,
		pattern='encode.*model',
		dat)

aic <- aic.table(models)
dir.tab(aic,'aic')

m <- models[[which.min(aic[,'aic'])]]

groups <- read.csv('out/groups.csv')

int <- read.csv('interactions.csv',row.names=1)[,3:4]
g <- graph_from_edgelist(as.matrix(int))
write.graph(g,format='dot',file='tmp.dot')
out <- mkdate('interactions','dot')
system2('sed',paste("'s/name/label/' tmp.dot >",out))
system2('dot',paste('-Tsvg -O',out))

int <- rbind(int,setNames(int[,2:1],names(int)))
int <- rbind(as.matrix(int),t(sapply(unique(groups$Condition),rep,2)))
int <- int[!duplicated(int),]
dir.csv(int,'interactions.csv')

int <- apply(int,1,paste,collapse='->')

# k <- get.k(3:53,m$dists,groups$Condition,int)
k <- get.k(3:53,m$dists,groups$Condition,int,'directed')
# k2 <- get.k(3:53,clusts$dists,groups$Condition,int,'directed')
dir.csv(k,'k')
k.max <- k[which.max(k[,'ES']),'k']

if(opts$all_k){
	sapply(3:53,function(k){
		dir.csv(res.unif(c(0.01,3),k,m$encoded,int,
				 groups$Condition,
				 opts$leiden_reps),
			paste0('leiden.k',as.character(k)))
		})
} else{
	leidens <- res.unif(c(0.01,3),k.max,m$encoded,
			    int,groups$Condition,
			    opts$leiden_reps)
	dir.csv(leidens,'leiden')
}
# k <- descend(c(3,53),gsea.prot,
#              m$dists,groups$Condition,int,
#              breaks=10,rate=1,maxiter=100,discrete=T,min=3)
# dir.csv(k$trace,'k')

#res <- res.unif(c(0.01,3),k.max,m$encoded,int,groups$Condition,10)
# 
# res <- get.res(c(0.01,3),10,1,100, 
#                k.max,m$encoded, 
#                int,groups$Condition,parallel=T)
# dir.csv(res$trace,'leiden'clustsclusts22)

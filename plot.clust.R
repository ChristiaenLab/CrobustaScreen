source("R/modelfns.R")
source("R/test.knn.R")
source("R/leiden.R")
source("R/clustplots.R")
source("R/plotfns.R")
source("R/io.R")
source("R/gene.network.R")

library(optparse)
library(igraph)
library(dirfns)
library(purrr)
library(umap)
library(ggpubr)

parser <- data.parser()
parser <- add_option(parser, c("-c", "--clust_dir"), 
		     action = "store",
		     default = Sys.Date(),
		     help = "Location of clustering output")
parse.env(parser)
list2env(read.clusts(clust_dir), globalenv())

umap.coords <- umap(encoded)$layout
colnames(umap.coords) <- c("UMAP1", "UMAP2")
dir.csv(umap.coords, "umap")

knn <- get.knn(dists, k, "plus")
enrichCond(groups$Condition,
	   as.matrix(as_adjacency_matrix(knn)),
	   "knn.network.fr", layout.fruchterman.reingold)
enrichCond(groups$Condition,
	   as.matrix(as_adjacency_matrix(knn)),
	   "knn.network")

hyper <- get.hyper(knn, groups$Condition)
g <- poisGraph(hyper)
networkPois(g, "hyper.k",
	    colfn = colorRamp2(c(0, max(E(g)$weight)), 
			       c("white", "red")))

statplot(leidens, paste0("leiden.k", as.character(k)))

dir.plot("knn")(plot.edge, umap.coords, knn)

#dir.f(ggexport)(dot.col("embedding2",
#			as.data.frame(encoded),
#			col = groups$Condition, "condition"),
#		filename = "embedding.pdf")

plots <- lapply(names(leidens)[2:7], dot.stat, leidens)

es <- dot.stat("ES", ks)
dir.f(ggexport)(ggarrange(plotlist = list(es), ncol = 3, nrow = 3),
		filename = "ES.pdf")
 
plots <- append(plots, list(es))
arrange.stats(plots, "optimization")

sel <- sapply(leidens[, c(2:6)], which.max)
sel["recall"] <- which(leidens[, 1] == 
		       max(leidens[leidens[, 4] == 
			   max(leidens[, 4]), 1]))

clusts <- clusts[, sel]
colnames(clusts) <- names(leidens[2:6])

dists <- as.matrix(as_adjacency_matrix(knn))
row.names(dists) <- as.numeric(1:nrow(dists))
colnames(dists) <- as.numeric(1:nrow(dists))

dir.plot('knn.clust')(plot.pt, umap.coords, knn, clusts[, 1])

clustcond <- function(cond, clusts = NULL, ...){
	dir.pdf(paste0('umap/', gsub("/", "_", cond)))
	plot.edge(umap.coords, knn, clusts)
	points(umap.coords[groups$Condition == cond, ], ...)
	dev.off()
}
sapply(unique(groups$Condition), clustcond, 
       clusts = clusts[, 1], 
       pch = 19, cex = 0.8, col = 1)


clustplots(encoded,
	   clusts[, 1],
	   groups$Phenotype,
	   dists, 'pheno')

clustplots(encoded,
	   clusts[, 1],
	   NULL,
	   dists, legend.ncol = 1, legend.cex = 1)

dir.f(quantHeatmap)(params, filename = "params",
		    split = clusts[, 1],
		    cell.w = 0.012, cell.h = 0.005,
		    show_row_names = F)
dir.f(quantHeatmap)(z, filename = "z",
		    split = clusts[, 1],
		    cell.w = 0.012, cell.h = 0.005,
		    show_row_names = F)
dir.f(quantHeatmap)(encoded, filename = "embedding",
		    split = clusts[, 1],
		    cell.w = 0.001, cell.h = 0.005,
		    show_row_names = F)

dir.f(clusthyper, 'out')(groups[, "Condition", drop = F], clusts[, 1])
dir.f(clusthyper, 'out')(as.data.frame(pheno), clusts[, 1],
			filename = 'pheno')

g <- gene.network(knn, resolution, groups$Condition, 
		  mode = 'directed')

edgelist <- cbind(as.data.frame(as_edgelist(g)), E(g)$weight)
edgelist <- do.call(rbind, lapply(split(edgelist, edgelist[, 1]), function(x) x[order(x[, 3], decreasing = T)[1:min(nrow(x), 5)], ]))
reduced <- graph_from_edgelist(as.matrix(edgelist[, -3]), F)
write.dot(reduced, 'gene_network')
graph.pdf('gene_network', reduced)

#sel <- E(g)$weight > quantile(E(g)$weight, 0.9)
#E(reduced)$weight <- E(g)$weight[sel]
#
#dir.f(networkPois, 'out')(reduced, 'gene_network',
#			 colfn = col.abs(E(reduced)$weight),
#			 title = 'modularity')
#
#dir.f(write_graph, 'file')(reduced, format = 'dot',
#			  filename = 'gene.network.dot')

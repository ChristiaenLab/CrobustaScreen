source("R/leiden.R")
source("R/heatmapfns.R")
source("R/gene.network.R")

source("R/hyper.R")
source("R/clust.params.R")

source("R/dirfns.R")
source("R/plotfns.R")
source("R/clustplots.R")
source("R/io.R")

library(ComplexHeatmap)
library(optparse)
library(igraph)
library(purrr)
library(umap)
library(ggpubr)

set.seed(42)

parser <- data.parser()
parser <- add_option(parser, c("-c", "--clust_dir"), 
		     action = "store",
		     #default = Sys.Date(),
		     default = "data",
		     help = "Location of clustering output")
parser <- add_option(parser, c("-s", "--clust_sel_method"), 
		     action = "store",
		     #default = "combined_score",
		     default = "ES",
		     help = "resolution selection parameter method for clusters. One of \"combined_score\", \"ES\", \"log2error\", \"mean_silhouette\", \"nclusts\".")
parse.env(parser)

# read data into global env
list2env(read.clusts(clust_dir), globalenv())

dir.pdf("param_cor", width=20, height=20)
draw(cor.hm(params, name="pearson correlation"))
dev.off()
dir.pdf("z_cor", width=20, height=20)
draw(cor.hm(z, name="pearson correlation"))
dev.off()
dir.pdf("embedding_cor", width=6, height=5)
draw(cor.hm(encoded, name="pearson correlation"))
dev.off()

dir.pdf("z_E_cor", width=8, height=20)
draw(xycor.hm(z, encoded, name="pearson correlation", cell.w=0.25))
dev.off()

row.names(encoded) <- row.names(params)
dir.csv(encoded, "embeddings")

umap.coords <- umap(encoded)$layout
row.names(umap.coords) <- row.names(params)
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
dir.f(ggexport)(ggarrange(plotlist = list(es), 
			  ncol = 3, nrow = 3),
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

condsel <- conds %in% c('Col9a1', 'Rab5', 'Rabep', 'Ddr', 'NDE', 'Rock')

dir.plot('knn.cond')(plot.pt, umap.coords, knn, 
					 conds, condsel, legendpos='bottomleft')
dir.plot('knn.clust')(plot.pt, umap.coords, knn, 
					  clusts[,clust_sel_method], legendpos='bottomleft')

clustcond <- function(cond, clusts = NULL, ...){
	dir.pdf(paste0('umap/edge/', gsub("/", "_", cond)))
	plot.edge(umap.coords, knn, clusts)
	points(umap.coords[groups$Condition == cond, ], ...)
	dev.off()

	dir.pdf(paste0('umap/point/', gsub("/", "_", cond)))
	plot.pt(umap.coords, knn, clusts)
	points(umap.coords[groups$Condition == cond, ], ...)
	dev.off()
}
sapply(unique(groups$Condition), clustcond, 
       clusts = clusts[,clust_sel_method], 
       pch = 1, cex = 0.8, col = 1)


clustplots(encoded,
	   clusts[,clust_sel_method],
	   groups$Phenotype,
	   dists, 'pheno')

clustplots(encoded,
	   clusts[,clust_sel_method],
	   NULL,
	   dists, legend.ncol = 1, legend.cex = 1)

save.hm <- function(mat, name, filename, cell.w = 0.120, cell.h=0.005, ...) {
	condfile <- sprintf('%s.cond', filename)
	clustfile <- sprintf('%s.clust', filename)

    # Height: rows * cell height + column labels (approx 0.1in/char) + padding
    max_cn <- max(nchar(colnames(mat)), 0)
    h <- nrow(mat) * cell.h + (max_cn * 0.1) + 4

    # Width Base: cols * cell width + padding
    w_base <- ncol(mat) * cell.w + 4

    # Cond Width: base + split labels
    split_cond <- groups$Condition
    max_split_cond <- max(nchar(as.character(unique(split_cond))), 0)
    w_cond <- w_base + (max_split_cond * 0.1)

    dir.pdf(condfile, path = "out", width = w_cond, height = h)
    draw(hm.cell(mat, name=name,
		    split = split_cond,
		    cell.w = cell.w, cell.h = cell.h,
		    show_row_names = F, 
		    row_title_rot = 0, ...))
    dev.off()

    # Clust Width: base + split labels
    split_clust <- clusts[,clust_sel_method]
    max_split_clust <- max(nchar(as.character(unique(split_clust))), 0)
    w_clust <- w_base + (max_split_clust * 0.1)

    dir.pdf(clustfile, path = "out", width = w_clust, height = h)
    draw(hm.cell(mat, name=name,
		    split = split_clust,
		    cell.w = cell.w, cell.h = cell.h,
		    show_row_names = F, 
		    row_title_rot = 0, ...))
    dev.off()
}

save.hm(params, "value", "params")
save.hm(z, "z-score",  "z")
save.hm(encoded, "value","embedding")

dir.f(clusthyper, 'out')(groups[, "Condition", drop = F], 
			 clusts[,clust_sel_method], 
			 filename = "condition")
dir.f(clusthyper, 'out')(as.data.frame(pheno), clusts[,clust_sel_method],
			filename = 'pheno')

dir.f(clustparam, "out")(params, clusts[,clust_sel_method],
			 filename = "params")
dir.f(clustparam, "out")(z, clusts[,clust_sel_method],
			 filename = "z")
dir.f(clustparam, "out")(encoded, clusts[,clust_sel_method],
			 filename = "embeddings")

dir.f(clustparam, "out")(params, groups$Condition,
			 filename = "condition/params")
dir.f(clustparam, "out")(z, groups$Condition,
			 filename = "condition/z")
dir.f(clustparam, "out")(encoded, groups$Condition,
			 filename = "condition/embeddings")

dir.f(clustparam, "out")(params, groups$Condition,
						 logfc.cutoff = 0.25, fdr.cutoff = 0.05, 
						 subset = c("Arhgef8", "Depdc", "Tyrosinase"),
						 filename = "Arhgef8_Depdc_Tyr/params")
dir.f(clustparam, "out")(z, groups$Condition,
						 logfc.cutoff = 0.25, fdr.cutoff = 0.05, 
						 subset = c("Arhgef8", "Depdc", "Tyrosinase"),
						 filename = "Arhgef8_Depdc_Tyr/z")
dir.f(clustparam, "out")(encoded, groups$Condition,
						 logfc.cutoff = 0, fdr.cutoff = 1, 
						 subset = c("Arhgef8", "Depdc", "Tyrosinase"),
						 filename = "Arhgef8_Depdc_Tyr/embeddings")

g <- gene.network(knn, resolution, groups$Condition, 
		  mode = 'directed')

edgelist <- cbind(as.data.frame(as_edgelist(g)), E(g)$weight)
edgelist <- do.call(rbind, 
    lapply(split(edgelist, 
		 edgelist[, 1]), 
	   function(x) { 
		   x[order(x[, 3], 
			   decreasing = T)[1:min(nrow(x), 5)], ]
	   }))

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


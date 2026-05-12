Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
options(reticulate.autoconfig = FALSE)
options(reticulate.conda_fallback = FALSE)

source("R/leiden.R")
source("R/gene.network.R")

source("R/hyper.R")
source("R/clust.params.R")

source("R/dirfns.R")
source("R/plotfns.R")
source("R/clustplots.R")
source("R/io.R")
source("R/anovas.R")

library(optparse)
library(igraph)
library(purrr)
library(umap)
library(ggpubr)

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

csvfn <- function(...) dir.csv(..., path = out_dir, append.date = F)
pdffn <- function(...) dir.pdf(..., path = out_dir, append.date = F)
writefn <- function(f, arg = "filename") {
	function(...) dir.f(f, arg)(..., path = out_dir, append.date = F)
}

omit <- c(-105:-106, -115:-117)
params <- params[, omit]

umap.coords <- umap(encoded)$layout
colnames(umap.coords) <- c("UMAP1", "UMAP2")

knn <- get.knn(dists, k, "plus")
dists <- as.matrix(as_adjacency_matrix(knn))
hyper <- get.hyper(knn, groups$Condition)
g <- poisGraph(hyper)

sel <- sapply(leidens[, c(2:6)], which.max)
sel["recall"] <- which(leidens[, 1] == 
               max(leidens[leidens[, 4] == 
               max(leidens[, 4]), 1]))

clusts <- clusts[, sel]
colnames(clusts) <- names(leidens[2:6])

csvfn(umap.coords, "umap")

enrichCond(groups$Condition,
       as.matrix(as_adjacency_matrix(knn)),
       paste0(out_dir, "/knn.network.fr"), layout.fruchterman.reingold)
enrichCond(groups$Condition,
       as.matrix(as_adjacency_matrix(knn)),
       paste0(out_dir, "/knn.network"))

networkPois(g, paste0(out_dir, "/hyper.k"),
        colfn = colorRamp2(c(0, max(E(g)$weight)), 
                   c("white", "red")))

statplot(leidens, paste0("leiden.k", as.character(k)), path = out_dir)

dir.plot("knn", pdffn)(plot.edge, umap.coords, knn)

#writefn(ggexport)(dot.col("embedding2",
#            as.data.frame(encoded),
#            col = groups$Condition, "condition"),
#        filename = "embedding.pdf")

plots <- lapply(names(leidens)[2:7], dot.stat, leidens)

if ("ES" %in% names(ks)) {
    es <- dot.stat("ES", ks)
    writefn(ggexport)(ggarrange(plotlist = list(es), 
                  ncol = 3, nrow = 3),
            filename = "ES.pdf")
    plots <- append(plots, list(es))
}
arrange.stats(plots, "optimization", path = out_dir)

write.aovs(params, z, encoded, groups, clusts[, clust_sel_method], out_dir)

row.names(dists) <- as.numeric(1:nrow(dists))
colnames(dists) <- as.numeric(1:nrow(dists))

dir.plot('knn.clust', pdffn)(plot.pt, umap.coords, knn, clusts[, clust_sel_method])

clustcond <- function(cond, clusts = NULL, ...) {
    pdffn(paste0('umap/edge/', gsub("/", "_", cond)))
    plot.edge(umap.coords, knn, clusts)
    points(umap.coords[groups$Condition == cond, ], ...)
    dev.off()

    pdffn(paste0('umap/point/', gsub("/", "_", cond)))
    plot.pt(umap.coords, knn, clusts)
    points(umap.coords[groups$Condition == cond, ], ...)
    dev.off()
}
sapply(unique(groups$Condition), clustcond, 
       clusts = clusts[, clust_sel_method], 
       pch = 1, cex = 0.8, col = 1)


clustplots(encoded,
       clusts[, clust_sel_method],
       groups$Phenotype,
       dists, 'pheno')

clustplots(encoded,
       clusts[, clust_sel_method],
       NULL,
       dists, legend.ncol = 1, legend.cex = 1)

writefn(quantHeatmap)(params[,omit], name = "value", filename = "params",
            split = clusts[, clust_sel_method],
            cell.w = 0.12, cell.h = 0.005,
            show_row_names = F)
writefn(quantHeatmap)(z, name = "z-score", filename = "z",
            split = clusts[, clust_sel_method],
            cell.w = 0.12, cell.h = 0.005,
            show_row_names = F)
writefn(quantHeatmap)(encoded, name = "value", filename = "embedding",
            split = clusts[, clust_sel_method],
            cell.w = 0.2, cell.h = 0.005,
            show_row_names = F)

writefn(clusthyper, 'path')(groups[, "Condition", drop = F],
              clusts[, clust_sel_method],
              filename = "clusters/condition")
writefn(clusthyper, 'path')(as.data.frame(pheno), clusts[, clust_sel_method],
             filename = 'clusters/pheno')
writefn(clustparam, "path")(params, clusts[, clust_sel_method],
             filename = "clusters/params", omit = omit)
writefn(clustparam, "path")(z, clusts[, clust_sel_method],
             filename = "z")
writefn(clustparam, "path")(encoded, clusts[, clust_sel_method],
             filename = "clusters/embeddings")

writefn(clusthyper, 'path')(as.data.frame(pheno), groups$Condition,
             filename = 'condition/pheno')
writefn(clustparam, "path")(params, groups$Condition,
              filename = "condition/params", omit = omit)
writefn(clustparam, "path")(z, groups$Condition,
              filename = "condition/z")
writefn(clustparam, "path")(encoded, groups$Condition,
              filename = "condition/embeddings")

writefn(clustparam, "path")(params, groups$Condition,
                         logfc.cutoff = 0.25, fdr.cutoff = 0.05,
                         subset = c("Arhgef8", "Depdc", "Tyrosinase"),
                         filename = "Arhgef8_Depdc_Tyr/params")
writefn(clustparam, "path")(z, groups$Condition,
                         logfc.cutoff = 0.25, fdr.cutoff = 0.05,
                         subset = c("Arhgef8", "Depdc", "Tyrosinase"),
                         filename = "Arhgef8_Depdc_Tyr/z")
writefn(clustparam, "path")(encoded, groups$Condition,
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
#writefn(networkPois, 'out')(reduced, 'gene_network',
#             colfn = col.abs(E(reduced)$weight),
#             title = 'modularity')
#
#writefn(write_graph, 'file')(reduced, format = 'dot',
#              filename = 'gene.network.dot')

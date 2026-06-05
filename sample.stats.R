# Compute per-sample silhouette width and within-cluster kNN connectivity
# for the optimal Leiden clustering selected by `--clust_sel_method`.

source("R/io.R")
source("R/dirfns.R")
source("R/sample.stats.R")

library(optparse)
library(igraph)

parser <- data.parser()
parser <- add_option(parser, c("-c", "--clust_dir"),
                     action = "store",
                     default = "data",
                     help = "Location of clustering output")
parser <- add_option(parser, c("-s", "--clust_sel_method"),
                     action = "store",
                     default = "ES",
                     help = paste("resolution selection method.",
                                  "One of \"combined_score\", \"ES\",",
                                  "\"recall\", \"log2error\", \"mean_silhouette\"."))
parse.env(parser)

list2env(read.clusts(clust_dir), globalenv())

# Pick the optimal clustering column. Mirrors the selection logic in
# plot.clust.R so the per-sample stats line up with the figures.
sel <- sapply(leidens[, c(2:6)], which.max)
sel["recall"] <- which(leidens[, 1] ==
                       max(leidens[leidens[, 4] ==
                           max(leidens[, 4]), 1]))

clust <- clusts[, sel[clust_sel_method]]

stats <- sample.stats(dists, clust, k,
                      sample.ids = groups$Row.names)

dir.csv(stats,
        paste0("sample_stats_", clust_sel_method),
        out_dir, append.date = FALSE)

exemplars <- cluster.exemplars(stats)

dir.csv(exemplars,
        paste0("cluster_exemplars_", clust_sel_method),
        out_dir, append.date = FALSE)

source("R/knn.R")

# Per-sample silhouette width and within-cluster kNN connectivity.
#
# `dists` is the full embedding-space distance matrix (as built by
# `read.embeddings`). `clust` is the integer cluster assignment vector,
# one entry per sample. `k` is the kNN size — should match the optimal
# `k` selected by `cluster.R` so that connectivity is measured on the
# same graph used everywhere else in the pipeline (`get.knn(..., "plus")`).
sample.stats <- function(dists, clust, k,
                         sample.ids = NULL,
                         knn.mode = "plus") {
    require(cluster)
    require(igraph)

    n <- length(clust)
    stopifnot(nrow(dists) == n, ncol(dists) == n)
    if (is.null(sample.ids)) sample.ids <- seq_len(n)

    # Per-sample silhouette width on embedding-space euclidean distance.
    sil <- silhouette(as.integer(clust), as.dist(dists))
    sil_width <- sil[, "sil_width"]

    # Within-cluster connectivity from the same kNN graph the rest of
    # the pipeline uses (cluster.R / plot.clust.R both call get.knn).
    g       <- get.knn(dists, k, knn.mode)
    adj_bin <- as.matrix(as_adjacency_matrix(g)) > 0
    diag(adj_bin) <- FALSE

    same <- outer(clust, clust, `==`)
    diag(same) <- FALSE

    deg   <- rowSums(adj_bin)
    intra <- rowSums(adj_bin & same)

    # Mean embedding-space distance from each sample to the rest of its
    # cluster. Singletons get NA. Lower = more centrally located in the
    # cluster.
    mean_intra_dist <- vapply(seq_len(n), function(i) {
        peers <- which(clust == clust[i])
        peers <- peers[peers != i]
        if (length(peers) == 0) NA_real_
        else mean(dists[i, peers])
    }, numeric(1))

    data.frame(sample          = sample.ids,
               cluster         = clust,
               sil_width       = sil_width,
               n_neighbors     = deg,
               n_intra_cluster = intra,
               frac_intra      = ifelse(deg > 0, intra / deg, NA_real_),
               mean_intra_dist = mean_intra_dist,
               row.names       = NULL)
}

# For each cluster, return the sample with the best value of each
# numeric metric in `stats` (the data frame produced by `sample.stats`).
# By default "best" means largest; metrics named in `minimize` are
# argmin'd instead (e.g. mean_intra_dist, where smaller = more central).
# Result has one row per cluster and one column per metric, with cells
# holding sample identifiers.
cluster.exemplars <- function(stats, metrics = NULL,
                              minimize = c("mean_intra_dist")) {
    if (is.null(metrics))
        metrics <- setdiff(names(stats), c("sample", "cluster"))

    clusters <- sort(unique(stats$cluster))
    out <- data.frame(cluster = clusters,
                      stringsAsFactors = FALSE)
    for (m in metrics) {
        pick <- if (m %in% minimize) which.min else which.max
        out[[m]] <- vapply(clusters, function(k) {
            ix <- which(stats$cluster == k)
            as.character(stats$sample[ix[pick(stats[[m]][ix])]])
        }, character(1))
    }
    out
}

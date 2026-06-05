#!/usr/bin/env bash

Rscript plot.clust.R -o "out/embeddings"
Rscript plot.clust.R -o "out/PCs" -e "data/PCs/E.csv" -c "data/PCs"
Rscript plot.clust.R -o "out/clusts1-3" -e "data/clusters_1-3/E.csv" -m "data/clusters_1-3" -c "data/clusters_1-3"

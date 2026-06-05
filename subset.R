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
parser <- add_option(parser, c("-n", "--n_clusts"),
	action="store",
	default="3", 
	help="number of clusters to write to `data/clusters1-{n_clusts}/E.csv`")
parse.env(parser)

# read data into global env
list2env(read.clusts(clust_dir), globalenv())

n <- as.numeric(n_clusts)
dir <- sprintf("data/clusters_1-%d", n)
clust.opts <- paste(" -o", dir, "-e", paste0(dir, "/E.csv"))
plot.opts <- paste(clust.opts, "-c", dir)

sel <- sapply(leidens[, c(2:6)], which.max)
sel["recall"] <- which(leidens[, 1] == 
		       max(leidens[leidens[, 4] == 
			   max(leidens[, 4]), 1]))
clusts <- clusts[, sel]
colnames(clusts) <- names(leidens[2:6])

sel <- clusts[,clust_sel_method] >= n

dir.csv(encoded[sel,], "E", dir, row.names=F, append.date=F)
dir.csv(groups[sel,], "groups", dir, row.names=F, append.date=F)
dir.csv(params[sel,], "params", dir, row.names=T, append.date=F)
dir.csv(z[sel,], "z_dat", dir, row.names=T, append.date=F)
dir.csv(pheno[sel,], "phenotype", dir, row.names=T, append.date=F)
file.copy(file.path(meta_dir, "interactions.csv"), 
	  file.path(dir, "interactions.csv"), overwrite = T)

clust.opts <- paste(" -o", dir, "-e", paste0(dir, "/E.csv"), "-m", dir)
plot.opts <- paste(clust.opts, "-c", dir)

system2("Rscript", paste("cluster.R", clust.opts))
system2("Rscript", paste("plot.clust.R", plot.opts))

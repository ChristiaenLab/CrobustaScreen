# search k and gamma values 
Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
options(reticulate.autoconfig = FALSE)
options(reticulate.conda_fallback = FALSE)

source('R/optimization.R')
source("R/io.R")
source("R/dirfns.R")

library(optparse)
library(igraph)
library(leiden)
#library(dirfns)

parser <- data.parser("data")

parser <- add_option(parser, c("-k", "--k_min"), 
     action = "store", 
     type = "integer",
     default = 3,
     help = "Minimum value of `k`")

parser <- add_option(parser, c("-K", "--k_max"), 
     action = "store", 
     type = "integer",
     default = 53,
     help = "Maximum value of `k`")

parser <- add_option(parser, c("-G", "--gamma_max"), 
     action = "store", 
     type = "double",
     default = 3.0,
     help = "Maximum value of `gamma`")

parser <- add_option(parser, c("-l", "--leiden_reps"), 
     action = "store", 
     type = "integer",
     default = 1000,
     help = "Number of times to run Leiden")

parser <- add_option(parser, c("-p", "--max_pcs"),
    action = "store", type = "integer", default = 50,
    help = "Maximum number of principal components to test")

parse.env(parser)

out <- paste0(out_dir, "/PCs")

pcs <- prcomp(z)$x

f <- function(n) {
    pcs <- pcs[,1:n]
	dists <- as.matrix(dist(pcs))

	get.k(k_min:k_max, dists, groups$Condition, interactions, 'directed')
	}

ks <- lapply(1:max_pcs, f)
ES <- do.call(cbind, lapply(ks,'[',,"ES"))
row.names(ES) <- ks[[1]][,"k"]
dir.csv(ES, 'k', out, append.date = F)

n_k = nrow(ES)
i <- which.max(ES)
i_r <- ((i - 1) %% n_k) + 1
i_c <- ((i - 1) %/% n_k) + 1

k <- (k_min:k_max)[i_r]
d <- i_c
E <- pcs[,1:d]
dir.csv(E, 'E', out, row.names = F, append.date = F)

leidens <- get.res.unif(c(0.05, gamma_max), k, E, 
    interactions, groups$Condition, 
    leiden_reps)

dir.csv(leidens, 'leiden', out, append.date = F)



system2("Rscript", paste("plot.clust.R", 
						 "-o", out, 
						 "-e", paste0(out, "/E.csv"),
						 "-c", out))

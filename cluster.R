# search k and gamma values 

source('R/optimization.R')
source("R/io.R")
source("R/dirfns.R")

library(optparse)
library(igraph)
#library(dirfns)

parser <- data.parser()

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

parse.env(parser)

ks <- get.k(k_min:k_max, dists, groups$Condition, interactions, 'directed')

dir.csv(ks, 'k', out_dir, append.date = F)

k <- ks[which.max(ks[,'ES']),'k']

leidens <- get.res.unif(c(0.01, gamma_max), k, encoded, 
			interactions, groups$Condition, 
			leiden_reps)

dir.csv(leidens, 'leiden', out_dir, append.date = F)

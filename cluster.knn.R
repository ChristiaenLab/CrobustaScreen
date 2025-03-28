# search k and gamma values 

source('R/optimization.R')
source("R/io.R")
source("R/dirfns.R")

library(optparse)
library(igraph)
library(leiden)
#library(dirfns)

parser <- data.parser("data")

parser <- add_option(parser, c("-k", "--k"), 
		     action = "store", 
		     type = "integer",
		     default = 3,
		     help = "value of `k`")

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

leidens <- get.res.unif(c(0.05, gamma_max), k, encoded, 
			interactions, groups$Condition, 
			leiden_reps)

dir.csv(leidens, 'leiden', out_dir, append.date = F)


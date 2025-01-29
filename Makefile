# Fetch putative orthologs between C. robusta, human, and mouse

all: data/interactions.csv data/X.csv

data/interactions.csv: STRINGdb
	Rscript get.interactions.R

STRINGdb: ensembl
	Rscript STRINGdb.R

ensembl: 
	Rscript cint.ensembl.R

data/X.csv: data/z_dat.csv
	julia preprocess.jl

data/params.csv:
	Rscript readPheno.R

data/embryodat.csv:
	Rscript readEmbryos.R

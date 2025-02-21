data.parser <- function(outdir = Sys.Date()){
	parser <- OptionParser()
	parser <- add_option(parser, 
			     c("-e", "--encoder_dir"), 
			     action = 'store',
			     default = 'data',
			     help = "Location of autoencoder output")
	parser <- add_option(parser, c("-m", '--meta_dir'), 
			     action = 'store',
			     default = 'data',
			     help = "Location of metadata")
	parser <- add_option(parser, c("-o", "--out_dir"), 
			     action = 'store',
			     default = outdir,
			     help = "Output directory")
	return(parser)
}

read.opt <- function(dir, file, ...){
	require(stringr)
	require(purrr)
	str_interp("${dir}/${file}") %>% 
		read.csv(...)
}

parse.interactions <- function(dir, groups){
int <- read.opt(dir, "interactions.csv",
		       row.names = 1)[,3:4]

	int <- rbind(int, setNames(int[,2:1], names(int)))
	int <- rbind(as.matrix(int), 
		     t(sapply(unique(groups$Condition),
			      rep,2)))
	int <- int[!duplicated(int),]

	int <- apply(int, 1, paste, collapse = '->')
	return(int)
}
	
read.params <- function(dir){
	require(purrr)

	groups <- read.opt(dir, "groups.csv")
	pheno <- read.opt(dir, "phenotype.csv",
			  row.names = 1)
	interactions <- parse.interactions(dir,groups)

	params <- read.opt(dir, "params.csv",
			  row.names = 1)
	z <- read.opt(dir, "z_dat.csv",
			  row.names = 1)

	colsel <- sapply(z,compose(abs,sum)) > 0
	params <- params[,colsel]
	z <- z[,colsel]

	list(groups = groups,
	     pheno = pheno,
	     params = params,
	     z = z,
	     interactions = interactions)
}

read.embeddings <- function(dir){
	encoded <- read.opt(dir, "E.csv")
	names(encoded) <- sub("Column", "embedding",
			      names(encoded))

	dists <- as.matrix(dist(encoded))
	list(encoded = encoded,
	     dists = dists)
}

read.clusts <- function(dir){
	ks <- read.opt(dir, "k.csv", row.names = 1)
	clusts <- read.opt(dir, "leiden.csv",
			   row.names = 1)
	clusts <- clusts[clusts$nclust > 1,]

	leidens <- clusts[,1:7]
	clusts <- t(clusts[,-1:-7])

	k <- ks[which.max(ks$ES),1]
	res <- leidens[which.max(leidens[,2]),1]
	list(clusts = clusts,
	     ks = ks,
	     leidens = leidens,
	     k = k,
	     resolution = res)
}

parse.env <- function(parser){
	opts <- parse_args(parser)
	
	dat <- read.params(opts$meta_dir)
	embedding <- read.embeddings(opts$encoder_dir)

	list2env(opts, globalenv())
	list2env(dat, globalenv())
	list2env(embedding, globalenv())
}

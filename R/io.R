data.parser <- function(outdir = Sys.Date()) {
  require(optparse)
  parser <- OptionParser()
  parser <- add_option(parser,
           #c("-e", "--encoder_dir"),
           c("-e", "--embeddings"),
           action = "store",
           default = "data/E.csv",
           help = "Location of autoencoder output")
  parser <- add_option(parser, c("-m", "--meta_dir"),
           action = "store",
           default = "data",
           help = "Location of metadata")
  parser <- add_option(parser, c("-o", "--out_dir"),
           action = "store",
           default = outdir,
           help = "Output directory")
  return(parser)
}

read.opt <- function(dir, file, ...) {
  require(stringr)
  require(purrr)
  str_interp("${dir}/${file}") %>%
    read.csv(...)
}

parse.interactions <- function(dir, groups) {
int <- read.opt(dir, "interactions.csv",
           row.names = 1)[,3:4]

  int <- rbind(int, setNames(int[,2:1], names(int)))
  int <- rbind(as.matrix(int),
         t(sapply(unique(groups$Condition),
            rep,2)))
  int <- int[!duplicated(int),]

  int <- apply(int, 1, paste, collapse = "->")
  return(int)
}
 
read.params <- function(dir) {
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

read.embeddings <- function(path, meta = NULL) {
  print(paste("Reading embeddings from:", path))
  encoded <- read.csv(path)

  # If the first column is character/factor, use it as row names
  if(!is.numeric(encoded[,1])) {
    rownames(encoded) <- encoded[,1]
    encoded <- encoded[,-1]
  }

  # Align with metadata if provided
  if(!is.null(meta)) {
    common_ids <- intersect(rownames(encoded), rownames(meta))
    if(length(common_ids) > 0) {
      encoded <- encoded[common_ids, , drop = FALSE]
    }
  }

  names(encoded) <- sub("Column", "embedding",
            names(encoded))

  # Ensure all data is numeric before dist()
  numeric_only <- encoded[sapply(encoded, is.numeric)]
  dists <- as.matrix(dist(numeric_only))
  
  list(encoded = numeric_only,
       dists = dists)
}

read.clusts <- function(dir){
	ks <- read.opt(dir, "k.csv", row.names = 1)
	clusts <- read.opt(dir, "leiden.csv",
			   row.names = 1)
	clusts <- clusts[clusts$nclust > 1,]

	leidens <- clusts[,1:7]
	clusts <- t(clusts[,-1:-7])

	# Handle ks as either a simple k,ES table or a matrix of scores
	if("ES" %in% names(ks)) {
		k <- ks$k[which.max(ks$ES)]
	} else {
		# If it's a matrix, find the max value across all columns
		best_idx <- which(ks == max(ks, na.rm = TRUE), arr.ind = TRUE)
		# which(..., arr.ind=TRUE) returns a matrix where rows are matches.
		# best_idx[1, 1] is the row index of the first match.
		k <- as.numeric(rownames(ks)[best_idx[1, 1]])
	}

	if(length(k) == 0 || is.na(k)) k <- as.numeric(rownames(ks)[1])

	res <- leidens[which.max(leidens[,2]),1]
	list(clusts = clusts,
	     ks = ks,
	     leidens = leidens,
	     k = k,
	     resolution = res)
}


parse.env <- function(parser) {
  opts <- parse_args(parser)

  dat <- read.params(opts$meta_dir)
  # Pass params to align embeddings
  embedding <- read.embeddings(opts$embeddings, meta = dat$params)

  list2env(opts, globalenv())
  list2env(dat, globalenv())
  list2env(embedding, globalenv())
}

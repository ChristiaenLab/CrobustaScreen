source("imageFns.R")
library(umap)
library(dirfns)

dirs <- list.files(
  "segdat",
  pattern = "_Statistics$",
  recursive = T,
  include.dirs = T,
  full.names = T
)

meta <- strsplit(dirs,'/')

meta <- data.frame(
  condition=sub('_','.',sub(' .*','',sapply(meta,'[',3))),
  grouping=dirs,
  batch=sapply(meta,'[',2),
  slide='',
  date=sub('.*([0-9]{2})([0-9]{2})([0-9]{4})','\\3_\\1_\\2',sapply(meta,'[',2)),
  image=paste0('image',sub('.*Image_([0-9]+)\\].*','\\1',dirs)),
  half='',
  file=dirs,
  stringsAsFactors = F
)
sel <- grep('slide',dirs)
meta$slide[sel] <- sub('.* (slide)\\s?([0-9]+) .*','\\1\\2',dirs[sel])

meta[grep('next part',dirs),'slide'] <- 'nextpart'
meta[grep('right',dirs,ignore.case=T),'half'] <- "right"
meta[grep('left',dirs,ignore.case=T),'half'] <- "left"

sel <- grep('normal',meta$batch)
meta[sel,'batch'] <- 'old'
meta[-sel,'batch'] <- 'new'
meta[meta$condition=="EPH",'condition'] <- "Eph"

meta$id <- sub('^_','',sub('_$','',paste(meta$slide,meta$image,meta$half,sep='_')))

meta[grep('grouped',meta[,2],ignore.case = T),2] <- 'surface'
meta[grep('surf[az]ce',meta[,2],ignore.case = T),2] <- 'surface'
meta[grep('individual',meta[,2],ignore.case = T),2] <- 'cell'
meta[grep('nuclei',meta[,2],ignore.case = T),2] <- 'cell'
meta[grep('cell',meta[,2],ignore.case = T),2] <- 'cell'

meta$name <- make.names(paste(meta$condition,meta$id,meta$date,sep='.'))

dat <- lapply(dirs,readEmbryos)

sel <- meta$grouping=='surface'
surface <- sel&!duplicated(meta$name[sel])
cell <- !sel&!duplicated(meta$name[!sel])
meta.surface <- meta[surface,]
meta.cell <- meta[cell,]
celldat <- dat[cell]
surfacedat <- dat[surface]

celldat <- lapply(celldat,normvol)

ndists <- sapply(celldat,getDists)
cdists <- sapply(celldat,getDists,F)
stats <- sapply(celldat,groupStats)
dists <- rbind(ndists,cdists,stats)

cellID <- sapply(celldat,function(x) factor(x[,2],levels=c("ATM","TVC")))
ncell <- sapply(cellID,table)
row.names(ncell) <- c('nATM','nTVC')
cell <- rbind(dists,ncell)
colnames(cell) <- meta.cell$name

membrane <- sapply(surfacedat,membraneStats)
colnames(membrane) <- meta.surface$name

params <- merge(t(cell),t(membrane),0)
params[sapply(params,function(x) !is.finite(x)&is.numeric(x))] <- 0
dir.csv(params, 'params', 'out', row.names=F, append.date=F)

z <- zdf(params[,-1])
z <- cbind(Row.names=params[,1],as.data.frame(z))
dir.csv(z, 'z_dat', 'out', append.date=F, row.names=F)

phenotype <- read.csv('phenotype.csv',stringsAsFactors=F)
names(phenotype)[1] <- 'Date'
phenotype$Date <- gsub('-','_',phenotype$Date)
phenotype$Condition <- gsub('\\s','',phenotype$Condition)
phenotype$Condition <- gsub('_','.',phenotype$Condition)
phenotype$Condition[phenotype$Condition=="EPH"] <- "Eph"
phenotype$ID.embryo <- tolower(phenotype$ID.embryo)
phenotype$ID.embryo <- sub("^[0-9]+_",'',phenotype$ID.embryo)
phenotype$Name <- make.names(paste(
	phenotype$Condition,
	phenotype$ID.embryo,
	phenotype$Date,
	sep='.'
))
# names(phenotype)[9] <- 'Phenotype'
phenotype$Phenotype <- ''
phenotype$migration.perturbed <- grepl('migration',phenotype$Comment.on.the.embryo)
phenotype$division.perturbed <- grepl('division',phenotype$Comment.on.the.embryo)
phenotype$orientation.perturbed <- grepl('orientation',phenotype$Comment.on.the.embryo)

phenotype$Phenotype <- sapply(
	phenotype$migration.perturbed,
	function(x) if(x) "migration" else ''
)
phenotype$Phenotype <- paste(phenotype$Phenotype,sapply(
	phenotype$division.perturbed,
	function(x) if(x) "division" else ''
),sep='/')
phenotype$Phenotype <- paste(phenotype$Phenotype,sapply(
	phenotype$orientation.perturbed,
	function(x) if(x) "orientation" else ''
),sep='/')
phenotype$Phenotype <- sub('^\\/+','',phenotype$Phenotype)
phenotype$Phenotype <- sub('\\/+$','',phenotype$Phenotype)
phenotype$Phenotype <- sub('\\/+','\\/',phenotype$Phenotype)
phenotype$Phenotype[phenotype$Phenotype!=''] <- paste(
	phenotype$Phenotype[phenotype$Phenotype!=''],
	'perturbed'
)
phenotype$Phenotype[phenotype$Phenotype==''] <- "WT"
phenotype <- phenotype[!duplicated(phenotype$Name),]

dat.phenotype <- merge(params,phenotype,by.x='Row.names',by.y='Name',all.x=T)
dat.phenotype$Condition <- sub('\\..*','',dat.phenotype[,1])
dat.phenotype$Date <- sub('.*\\.','',dat.phenotype[,1])

dir.csv(dat.phenotype,'dat','out',row.names=F,append.date=F)

groups <- dat.phenotype[,c(
	'Row.names', 'Condition', 'Phenotype', #'nTVC', 'nATM',
	'migration.perturbed', 'orientation.perturbed', 'division.perturbed'
	#'TVC.contiguous', 'ATM.contiguous'
)]
dir.csv(groups,'groups','out',row.names=F,append.date=F)

pheno.sub <- groups[
	groups$nTVC==4&groups$nATM==2&groups$TVC.contiguous&groups$ATM.contiguous,
]
dir.csv(pheno.sub,'phenoSub','out',row.names=F,append.date=F)
dir.csv(params[params[,1]%in%pheno.sub[,1],],'phenoSubParam','out',row.names=F,append.date=F)
dir.csv(z[z[,1]%in%pheno.sub[,1],],'phenoSubZ','out',row.names=F,append.date=F)

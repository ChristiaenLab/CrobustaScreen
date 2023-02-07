library(dirfns)
dat <- read.csv('out/z_dat.csv',row.names=1)

phenotype <- read.csv('phenotype.csv',stringsAsFactors=F)
names(phenotype)[1] <- 'Date'
phenotype$Date <- gsub('-','_',phenotype$Date)

phenotype$Condition <- sub('.*\\s','',phenotype$Condition)
phenotype$Condition[phenotype$Condition=="EPH"] <- "Eph"
phenotype$Condition[phenotype$Condition=="GNA-L/S"] <- "Gna-L/S"

phenotype$Condition[phenotype$Condition=="12/13"] <- "Gna12/13"
phenotype$ID.embryo <- tolower(phenotype$ID.embryo)
phenotype$ID.embryo <- sub("^[0-9]+_",'',phenotype$ID.embryo)

phenotype$ID.embryo <- sub("with_lightning_","",phenotype$ID.embryo)

phenotype$Name <- make.names(paste(
	phenotype$Condition,
	phenotype$ID.embryo,
	phenotype$Date,
	sep='.'
))

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
row.names(phenotype) <- make.unique(phenotype$Name)


tmp <- setdiff(params$Row.names,row.names(phenotype))
tmp2 <- setdiff(row.names(phenotype),params$Row.names)

sel <- sapply(sub('\\.[0-9]{4}_{0-9}{2}_{0-9}{2}.*','',tmp),grep,row.names(phenotype),ignore.case=T)
row.names(phenotype)[sel] <- tmp

dat.phenotype <- merge(params,phenotype,by.x='Row.names',by.y=0,all.x=T)
dat.phenotype$Date <- sub('.*\\.','',dat.phenotype[,1])

dir.csv(dat.phenotype,'dat','out',row.names=F,append.date=F)

groups <- dat.phenotype[,c(1,(length(params)+1):length(dat.phenotype))]

pheno.sub <- groups[
	dat.phenotype$nTVC==4&dat.phenotype$nATM==2&dat.phenotype$TVC.contiguous&dat.phenotype$ATM.contiguous,
]

dir.csv(groups,'groups', 'out', row.names=F, append.date=F)
d

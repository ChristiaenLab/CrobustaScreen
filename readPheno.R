params <- read.csv("data/embryodat.csv")

phenotype <- read.csv('data/imaris.csv', stringsAsFactors = F)
names(phenotype)[1] <- 'Date'
phenotype$Date <- gsub('-','_',phenotype$Date)
# phenotype$Condition <- gsub('\\s','',phenotype$Condition)
# 
# phenotype$Cond <- sub('^[0-9]+[A-Z][-\\/]','',phenotype$Cond)
# phenotype$Cond <- sub("Gna12\\/[0-9]+","Gna12.13",phenotype$Cond)
# phenotype$Cond <- sub("EPB4\\.[0-9]+","EPB4.1",phenotype$Cond)
# 
# phenotype$Cond <- gsub('_','.',phenotype$Cond)

phenotype$Condition <- sub('.*\\s','',phenotype$Condition)
phenotype$Condition[phenotype$Condition=="EPH"] <- "Eph"
phenotype$Condition[phenotype$Condition=="GNA-L/S"] <- "Gna-L/S"
# phenotype$Condition[phenotype$Condition=="GNA12.13"] <- "Gna12.13"
phenotype$Condition[phenotype$Condition=="12/13"] <- "Gna12/13"
phenotype$ID.embryo <- tolower(phenotype$ID.embryo)
phenotype$ID.embryo <- sub("^[0-9-]+_",'',phenotype$ID.embryo)
#phenotype <- phenotype[-grep("with_lightning",phenotype$ID.embryo)]
phenotype$ID.embryo <- sub("with_lightning_","",phenotype$ID.embryo)
phenotype$Threshold.membrane <- as.numeric(phenotype$Threshold.membrane)

phenotype$Name <- make.names(paste(
	phenotype$Condition,
	phenotype$ID.embryo,
	phenotype$Date,
	sep='.'
))

# names(phenotype)[9] <- 'Phenotype'
phenotype$phenotype.label <- ''
phenotype$migration.perturbed <- grepl('migration',phenotype$Comment.on.the.embryo)
phenotype$division.perturbed <- grepl('division',phenotype$Comment.on.the.embryo)
phenotype$orientation.perturbed <- grepl('orientation',phenotype$Comment.on.the.embryo)

phenotype$phenotype.label <- sapply(
	phenotype$migration.perturbed,
	function(x) if(x) "migration" else ''
)
phenotype$phenotype.label <- paste(phenotype$phenotype.label,sapply(
	phenotype$division.perturbed,
	function(x) if(x) "division" else ''
),sep='/')
phenotype$phenotype.label <- paste(phenotype$phenotype.label,sapply(
	phenotype$orientation.perturbed,
	function(x) if(x) "orientation" else ''
),sep='/')
phenotype$phenotype.label <- sub('^\\/+','',phenotype$phenotype.label)
phenotype$phenotype.label <- sub('\\/+$','',phenotype$phenotype.label)
phenotype$phenotype.label <- sub('\\/+','\\/',phenotype$phenotype.label)
phenotype$phenotype.label[phenotype$phenotype.label!=''] <- paste(
	phenotype$phenotype.label[phenotype$phenotype.label!=''],
	'perturbed'
)
phenotype$phenotype.label[phenotype$phenotype.label==''] <- "WT"
row.names(phenotype) <- make.unique(phenotype$Name)

#phenotype <- phenotype[!duplicated(phenotype$Name),]

tmp <- setdiff(params$Row.names,row.names(phenotype))
tmp2 <- setdiff(row.names(phenotype),params$Row.names)

sel <- sapply(sub('\\.[0-9]{4}_{0-9}{2}_{0-9}{2}.*','',tmp),grep,row.names(phenotype),ignore.case=T)
ambiguous <- sapply(sel,length) > 1

row.names(phenotype)[unlist(sel[!ambiguous])] <- tmp[!ambiguous]
sel <- sel[ambiguous]
tmp <- tmp[ambiguous]

M <- lapply(sel,function(i) phenotype[i,c("migration.perturbed","division.perturbed","orientation.perturbed")])
same <- sapply(M,function(m) all(sapply(m,function(x) length(unique(x))==1)))

row.names(phenotype)[sapply(sel[same],`[`,1)] <- tmp[same]

tmp <- tmp[!same]
sel <- sel[!same]
M <- M[!same]

wt <- params[params$Row.names %in% tmp, c("nATM","nTVC")]
wt <- wt$nATM == 2 & wt$nTVC == 4
wt.y <- lapply(M, function(x) x[,1] & x[,2])
cells.match <- mapply(`==`, wt, wt.y, SIMPLIFY=F)
ambiguous <- sapply(cells.match,sum) != 1
i <- sapply(cells.match[!ambiguous],which)
j <- mapply(`[[`,sel[!ambiguous],i,SIMPLIFY=T)

row.names(phenotype)[j] <- tmp[!ambiguous]

dat.phenotype <- merge(params,phenotype,by.x='Row.names',by.y=0,all.x=T)
# dat.phenotype$Condition <- sub('\\..*','',dat.phenotype[,1])
dat.phenotype$Date <- sub('.*\\.','',dat.phenotype[,1])

dir.csv(dat.phenotype,'dat','out',row.names=F,append.date=F)

groups <- dat.phenotype[,c(1,(length(params)+1):length(dat.phenotype))]
			   #         'Row.names', 'Condition', 'Phenotype',# 'nTVC', 'nATM',
			   #         'migration.perturbed', 'orientation.perturbed', 'division.perturbed'
	#         'TVC.contiguous', 'ATM.contiguous'
# )]
dir.csv(groups,'groups', 'out', row.names=F, append.date=F)

pheno.sel <- sapply(compose(unique, unlist, 
			strsplit)(groups$Phenotype,
			split=' / '),
		grep,
		groups$Phenotype)

pheno <- replicate(5,rep('WT',nrow(groups)))
colnames(pheno) <- c("TVC.division","TVC.migration",
		     "ATM.division","ATM.migration",
		     "other")
row.names(pheno) <- row.names(groups)
pheno[pheno.sel$`inhibited TVC division`,
      "TVC.division"] <- "inhibited"
pheno[pheno.sel$`enhanced TVC division`,
      "TVC.division"] <- "enhanced"
pheno[pheno.sel$`inhibited ATM division`,
      "ATM.division"] <- "inhibited"
pheno[pheno.sel$`enhanced ATM division`,
      "ATM.division"] <- "enhanced"
pheno[pheno.sel$`inhibited TVC migration`,
      "TVC.migration"] <- "inhibited"
pheno[pheno.sel$`enhanced TVC migration`,
      "TVC.migration"] <- "enhanced"
pheno[pheno.sel$`inhibited ATM migration`,
      "ATM.migration"] <- "inhibited"
pheno[pheno.sel$`enhanced ATM migration`,
      "ATM.migration"] <- "enhanced"
pheno[pheno.sel$`problem disposition TVC`,
      "other"] <- "TVC disposition"
pheno[pheno.sel$`TVC cell alignment`,
      "other"] <- "TVC alignment"
# pheno[pheno.sel$`problem disposition ATM`,
#       "other"] <- "ATM disposition"
# pheno <- cbind(Condition=groups$Condition,pheno)
row.names(pheno) <- groups$Row.names
dir.csv(pheno,'phenotype', 'out', append.date=F)


pheno.sub <- groups[
	dat.phenotype$nTVC==4&dat.phenotype$nATM==2&dat.phenotype$TVC.contiguous&dat.phenotype$ATM.contiguous,
]

M <- sapply(dat.phenotype[,c("Threshold.nuclei","Threshold.membrane",
		    "Shortest.distance.between.ATM.group.to.TVC.group..from.edge.to.edge..µm")],as.numeric)
colnames(M)[3] <- "surface.dist"
params <- cbind(params,as.data.frame(M))

params[sapply(params,function(x) !is.finite(x)&is.numeric(x))] <- 0
dir.csv(params, 'params', 'out', row.names=F, append.date=F)

z <- zdf(params[,-1])
z <- cbind(Row.names=params[,1],as.data.frame(z))
dir.csv(z, 'z_dat', 'out', append.date=F, row.names=F)

dir.csv(groups,'groups', 'out', row.names=F, append.date=F)
dir.csv(pheno.sub,'phenoSub','out',row.names=F,append.date=F)
dir.csv(params[params[,1]%in%pheno.sub[,1],],'phenoSubParam','out',row.names=F,append.date=F)
dir.csv(z[z[,1]%in%pheno.sub[,1],],'phenoSubZ','out',row.names=F,append.date=F)



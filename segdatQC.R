library(moreComplexHeatmap)
source("imageFns.R")

dirs <- list.files(
  "segdat_adj",
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
  date=sub('.*([0-9]{2})([0-9]{2})([0-9]{4})_.*','\\3_\\1_\\2',sapply(meta,'[',3)),
  image=paste0('image',sub('.*Image_([0-9]+)\\].*','\\1',dirs)),
  half='',
  file=dirs,
  stringsAsFactors = F
)

sel <- grep('slide',dirs)
meta$slide[sel] <- sub('.* (slide)\\s?([0-9]+) .*','\\1\\2',dirs[sel])

meta[grep('right',dirs,ignore.case=T),'half'] <- "right"
meta[grep('left',dirs,ignore.case=T),'half'] <- "left"

meta[meta$condition=="depdc",'condition'] <- "Depdc"

meta$id <- sub('^_','',sub('_$','',paste(meta$slide,meta$image,meta$half,sep='_')))

meta$name <- make.names(paste(meta$condition,meta$id,meta$date,sep='.'))

dat <- lapply(dirs,readEmbryos)

dat <- split(dat,meta$batch)

meta <- split(meta,meta$batch)

parsed <- lapply(dat,function(x) sapply(x, membraneStats)[-7:-8,])

mapply(
	cor, 
	as.data.frame(t(parsed$adj)), 
	as.data.frame(t(parsed$unadj)), 
	use='na.or.complete'
)


dir.eps('qc')
hm <- Heatmap(
	t(rbind(parsed$adj,parsed$unadj)),
	column_split=c(rep('Adjusted',nrow(parsed$adj)), rep('Unadjusted',nrow(parsed$unadj))),
	right_annotation=rowAnnotation(condition=meta$adj$condition)
)
draw(hm)
dev.off()

mapply(wilcox.test,as.data.frame(parsed[[1]]),as.data.frame(parsed[[2]]))
mapply(t.test,as.data.frame(t(parsed$adj)),as.data.frame(t(parsed$unadj)))

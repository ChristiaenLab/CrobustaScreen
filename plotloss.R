source("R/heatmapfns.R")

library(ComplexHeatmap)
library(purrr)

loss <- read.csv("data/loss_DEWAK.csv")
loss <- loss[!duplicated(loss[,1:2]),]

dat <- reshape(loss, v.names = names(loss)[-1:-2], 
			   timevar = "d", idvar = "k", direction = "wide")
row.names(dat) <- dat$k
dat <- dat[,-1]
sp  <- sub("\\..*", "", names(dat))

Ms <- lapply(split(as.data.frame(t(dat)), sp), 
			 function(M){ 
				 row.names(M) <- sub(".*\\.", "", row.names(M))
				 t(M)
			 })
names(Ms) <- unique(sp)

hms <- mapply(function(M, name, col) {
				  Heatmap(M, col.abs(M, cols = c("white",col)), name, 
						  row_title = "k", column_title = "d")
			 }, Ms, names(Ms), c("black", "blue", "red", "darkgreen"))
hm <- Reduce(`+`, hms)

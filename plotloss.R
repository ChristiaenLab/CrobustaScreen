source("R/heatmapfns.R")
source("R/plotfns.R")

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
labs <- unique(sp)


hms <- mapply(function(M, name, col) {
				  Heatmap(M, col.abs(M, cols = c("white",col)), name, 
						  row_title = "k", column_title = "d",
				  		  cluster_rows = F, cluster_columns = F)
			 }, Ms, labs, c("black", "blue", "red", "darkgreen"))
hm <- Reduce(`+`, hms)

dir.pdf("loss_DEWAK.pdf")
draw(hm)
dev.off()

loss <- read.csv("data/loss_dk.csv")
loss_d <- loss[1:54, -2]
loss_k <- loss[55:nrow(loss), -1]
plots_d <- lapply(names(loss)[3:6], dot.stat, loss_d)
plots_k <- lapply(names(loss)[3:6], dot.stat, loss_k)
plots <- append(plots_d, plots_k)
arrange.stats(plots, "loss_dk", ncols = 4, nrows = 4)

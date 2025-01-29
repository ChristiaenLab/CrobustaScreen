source("R/hyper.R")

silhouettePlot <- function(dat, clusts, distmat = NULL, out){
	require(dirfns)
	require(cluster)
	pcs <- dat[, grep('pca', names(pcs))]
	pcs <- dat[, !apply(pcs == 0, 2, all)]

	row.names(pcs) <- row.names(dat)

	if(is.null(distmat)) {
		   dists <- dist(pcs)
	}else dists <- distmat

	sil <- silhouette(clusts, dists)

#	if(!is.na(sil)){
		dir.pdf('silhouette', out, append.date = F)
		plot(sil)
		dev.off()
#	}
	return(sil)
}

plotDist <- function(dat, conds, distmat, out, clustcol = 'gray', ...){
	require(dirfns)

	dat <- dat[, grep('umap', names(dat))]
	cols <- c(Depdc = 'royalblue', Tyrosinase = 'saddlebrown')
	condcol <- cols[conds]
	#         dists <- distmat[row.names(dat), row.names(dat)]

	ix <- which(distmat > 0, arr.ind = T)
	clustline <- do.call(rbind, 
			     apply(ix, 1, 
				   function(x) rbind(dat[x[1], ], dat[x[2], ])))

	dir.pdf('Tyr_Depdc', out, append.date = F)
	plot(NULL,
	     xlim = range(dat[, 1]), ylim = range(dat[, 2]), 
	     xlab = NA, ylab = NA, xaxt = 'n', yaxt = 'n')
	lines(clustline, col = clustcol)
	points(dat, cex = 1.0, col = condcol, pch = 19)
	legend(
	       'bottomleft',
	       names(cols), 
	       col = cols,
	       pch = 19,
	       cex = 1.0
	)
	dev.off()
}

clustplots <- function(dat, clusts, conds, distmat, 
		       out = '.', append.date = T,
		       legend.ncol = 2, legend.cex = 0.5){
	require(dirfns)
	require(cluster)

	row.names(dat) <- as.character(1:nrow(dat))
	umap <- as.data.frame(dat)
	if(length(clusts) > 0) {
		clustpts <- split(umap, clusts)
	}else clustpts <- NULL
	if(length(conds) > 0) {
		condpts <- split(umap, conds)
	} else condpts <- NULL

	clustcol <- rainbow(length(clustpts))
	condcol <- rainbow(length(condpts))

	clustdist <- lapply(clustpts, 
			    function(x) {
				    distmat[row.names(x), 
					    row.names(x)]
			    })

	clustline <- mapply(function(dists, pts){
		ix <- which(dists > 0, arr.ind = T)
		do.call(rbind, 
			apply(ix, 1, 
			      function(x) { 
				      rbind(pts[x[1], ], 
					    pts[x[2], ])
			      }))
	}, clustdist, clustpts, SIMPLIFY = F)

	dir.pdf('clustedge', out, append.date = append.date)
	plot(NULL, xlim = range(umap[, 1]), 
	     ylim = range(umap[, 2]), xlab = NA, ylab = NA, 
	     xaxt = 'n', yaxt = 'n')
	mapply(lines, clustline, col = clustcol)
	mapply(points, condpts, col = condcol, cex = 0.5)
	legend(
	       'bottomleft',
	       c(names(condpts), names(clustpts)),
	       col = c(condcol, clustcol),
	       lty = c(rep(NA, length(condpts)), 
		     rep(1, length(clustpts))),
	       pch = c(rep(1, length(condpts)), 
		     rep(NA, length(clustpts))),
	       ncol = legend.ncol, cex = legend.cex
	)
	dev.off()

	dir.pdf('clustplot', out, append.date = append.date)
	plot(NULL, xlim = range(umap[, 1]), 
	     ylim = range(umap[, 2]), xlab = NA, ylab = NA, 
	     xaxt = 'n', yaxt = 'n')
	mapply(points, condpts, col = condcol, cex = 0.5, pch = 19)
	mapply(points, clustpts, col = clustcol, cex = 0.5)
	legend(
	       'bottomleft',
	       c(names(condpts), names(clustpts)),
	       col = c(condcol, clustcol),
	       pch = c(rep(19, length(condpts)), 
		     rep(1, length(clustpts))),
	       ncol = 2, cex = 0.5
	)
	dev.off()
}

clustparam <- function(dat, clusts, out){
	require(moreComplexHeatmap)
	clustdat <- split(dat, clusts)
	# test for significance of each feature in each cluster
	testfn <-function(x, clust, dat){
	      mu <- mean(dat[, x], na.rm = T)
	      FC <- mean(clust[, x], na.rm = T)/mu
	      if(FC != 1 & length(unique(clust[, x])) > 1){
		      utest <- wilcox.test(clust[, x], mu = mean(dat[, x], na.rm = T))$p.value
		      ttest <- t.test(clust[, x], mu = mean(dat[, x], na.rm = T))$p.value
	      } else {
		      utest <- NaN
		      ttest <- NaN
	      }
	      return(c(FC = FC, u = utest, t = ttest))
	} 
	clusttest <- function(clust){
	       sapply(colnames(clust), testfn, clust, dat)
	}
	test <- lapply(clustdat, clusttest)

	# select fields from output
	fcdat <- sapply(test, '[', "FC", )
	p.t <- sapply(test, '[', 't', )
	p.u <- sapply(test, '[', 'u', )

	# convert p-values to FDR values
	fdr.t <- apply(p.t, 2, function(x) p.adjust(unlist(x)))

	# define color scale
	clustdat <- lapply(clustdat, function(x) t(x)[row.names(fcdat), ])
	# col.z creates a color scale centered on 0 ranging from the 0.01 to 0.99 quantiles of the input data
	if(length(unique(clusts)) > 1){
	# This scale will be used for the fold changes between average feature values in each cluster from the background.
		colfc <- col.z(fcdat, .05, 1)
		# creates a color scale from 0 to 2
		# This will give a log10 scale for FDR values between 1 and 0.01
		colfdr <- colorRamp2(c(0, 2), c('white', 'black'))

		# show feature values within clusters as boxplots color-coded by FDR value
		getAnn <- function(x) {
			m <- clustdat[[x]]
			fc.cols <- colfc(fcdat[, x])
			fdr.cols <- colfdr(fdr.t[, x])
			return(anno_boxplot(
				m, 
				which = 'row',
				width = unit(1, 'in'),
				box_width = 0.9,
				gp = gpar(
					fill = fc.cols,
					col = fdr.cols
				)
			))
		}

		# apply annotation function to each cluster
		ha <- lapply(colnames(fcdat), getAnn)
		names(ha) <- paste('cluster', colnames(fcdat))
		# bind annotation into single object
		ha <- do.call(HeatmapAnnotation, append(ha, list(which = 'row')))

		# assign names to color scales
		lgd <- list(
			Legend(col_fun = colfc, title = "FC"),
			Legend(col_fun = colfdr, title = "-log10(FDR)")
		)

		hm <- hm.cell(fcdat,
			right_annotation = ha,
			cell.w = .15,
			cell.h = .15,
			show_heatmap_legend = F,
			col = colfc
		)

		dir.pdf('t', out, height = 24, width = 10+ncol(fcdat), append.date = F)
		draw(hm, annotation_legend_list = lgd)
		dev.off()
	}
}

clusthyper <- function(dat, clusts, out){
	require(moreComplexHeatmap)
	# run hypergeometric tests for enrichment of conditions and phenotypes
	hyper <- lapply(dat, 
			function(x) condHyper(row.names(dat),
					      x, clusts))

	# extract fields from test & reformat as matrices
	odds <- do.call(rbind, lapply(hyper, '[[', 'log2OR'))
	fdr <- do.call(rbind, lapply(hyper, '[[', 'FDR'))
	qval <- as.matrix(do.call(rbind,
				  lapply(hyper, '[[', 'q')))

	# split phenotype & condition into separate panels
	#         if(length(dat) > 1){
	rowsplit <- unlist(
		mapply(
		       function(x, y) rep(y, nrow(x$log2OR)),
		       hyper,
		       names(hyper)
		)
	)
	#         } else rowsplit <- NULL

	# write matrices to csv
	dir.csv(cbind(
		parameter = rowsplit,
		condition = row.names(odds),
		as.data.frame(odds)
	), 'log2OR', out, append.date = F)
	dir.csv(cbind(
		parameter = rowsplit,
		condition = row.names(odds),
		as.data.frame(fdr)
	), 'FDR', out, append.date = F)
	dir.csv(cbind(
		parameter = rowsplit,
		condition = row.names(odds),
		as.data.frame(qval)
	), 'size', out, append.date = F)

	if(length(unique(clusts)) > 1){
		dotPscale(
			odds, 
			fdr, 
			qval, 
			file = 'condition', 
			path = out, 
			row_split = rowsplit, 
			row_title_rot = 0
		)
	}
}

condDist <- function(clust, cond, dists, out){
	require(moreComplexHeatmap)

	condsel <- sapply(unique(cond), function(x) cond == x)
	conds <- combn(unique(cond), 2)
	conds <- cbind(sapply(unique(cond), rep, 2), conds)
	sel <- sapply(unique(clust), function(x) clust == x)
	colnames(sel) <- unique(clust)

	tmp <- lapply(as.data.frame(sel), function(x){
		x <- setNames(apply(conds, 2, function(y){
			as.data.frame(dists[x & condsel[, y[1]],
				     x & condsel[, y[2]], drop = F])
		}), apply(conds, 2, paste, collapse = '.'))
		x <- lapply(x, function(x) x[x != 0])
		x
	})

	mat <- sapply(tmp, sapply, mean, na.rm = T)
	sp <- sub('\\..*', '', row.names(mat))

	avg <- apply(conds, 2, function(x) mean(dists[condsel[, x[1]], condsel[, x[2]]]))
	avg <- split(avg, sp)
	avg <- do.call(c, sapply(avg, '[', -1))

	ct <- sapply(tmp, sapply, length)
	ct <- split(as.data.frame(ct), sp)
	ct <- do.call(rbind, lapply(ct, '[', -1, ,drop = F))

	dat <- split(as.data.frame(mat), sp)
	lapply(dat, function(x) as.matrix(x)/x[1, ])
	self <- lapply(tmp, function(x) x[1:length(unique(cond))])
	self <- sapply(self, sapply, mean, na.rm = T)
	row.names(self) <- sub('\\..*', '', row.names(self))
	dat <- sapply(names(dat), function(x) t(t(dat[[x]])/self[x, ])[-1, ,drop = F])
	dat <- do.call(rbind, dat)
	sp <- sub('\\..*', '', row.names(dat))
	row.names(dat) <- sub('.*\\.', '', row.names(dat))

	n <- table(cond)
	n <- paste0(names(n), ' (', as.character(n), ')')
	sp <- Reduce(function(x, y) sub(sub('\\s.*', '', y), y, x), n, sp)

	size <- apply(sel, 2, sum)
	score <- apply(dat*ct, 2, sum, na.rm = T)/apply(ct, 2, sum)
	total <- log2(sum(score*size)/sum(size))

	#         topann <- list(anno_barplot(size), anno_barplot(score))
	#         names(topann) <- c('size', paste('score =', as.character(total)))
	#         topann <- do.call(columnAnnotation, topann)

	topann <- columnAnnotation(size = anno_barplot(size),
				   score = anno_barplot(score),
				   name = paste("score =", as.character(total)))
	quantHeatmap(
		log2(dat), 'cond_dist', path = out,
		conds = row.names(dat), show_row_names = F,
		top_annotation = topann,
		right_annotation = rowAnnotation(avg.dist = anno_barplot(avg)),
		column_title = paste('total score =', as.character(total)),
		split = sp, row_title_rot = 0
	)

	out <- lapply(tmp, function(x) x[length(unique(cond))+1:length(x)])
	mean(unlist(self))/mean(unlist(out))
}

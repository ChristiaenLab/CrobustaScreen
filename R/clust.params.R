source("R/dotplot.R")

clust.t.test <- function(x, clust, dat){
        mu <- mean(dat[, x], na.rm = T)
        mu_cl <- mean(clust[, x], na.rm = T)
        FC <- mu_cl / mu
        if(FC != 1 & length(unique(clust[, x])) > 1){
	utest <- wilcox.test(clust[, x], 
			     mu = mu)$p.value
	ttest <- t.test(clust[, x], 
			mu = mu)$p.value
        } else {
                utest <- NaN
                ttest <- NaN
        }
        return(c(mu = mu_cl, FC = FC, 
		 u = utest, t = ttest))
} 

test.clust.params <- function(dat,clustdat){
# test for significance of each feature in each cluster
	clusttest <- function(clust){
	       sapply(colnames(clust), 
		      clust.t.test, 
		      clust, dat)
	}
	test <- lapply(clustdat, clusttest)
        return(test)
}

box.heatmap <- function(m, clustdat, boxdat, outldat,
		  out, path, 
		  boxtitle = "log2(FC)",
		  outltitle = "-log10(FDR)", ...){
	# This scale will be used for the fold changes 
	# between average feature values in each cluster 
	# from the background.
        boxcol <- col.z(boxdat, .05, 1)

        # creates a color scale from 0 to 2
        # This will give a log10 scale for FDR values 
	# between 1 and 0.01
        outlcol <- colorRamp2(c(0, 2), c('white', 'black'))

        # show feature values within clusters as boxplots 
	# color-coded by FDR value
        getAnn <- function(x) {
                m <- clustdat[[x]]
                fc.cols <- boxcol(boxdat[, x])
                fdr.cols <- outlcol(outldat[, x])
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
        ha <- lapply(colnames(boxdat), getAnn)
        names(ha) <- paste('cluster', colnames(boxdat))

        # bind annotation into single object
        ha <- do.call(HeatmapAnnotation, 
		      append(ha, list(which = 'row')))

        # assign names to color scales
        lgd <- list(
                Legend(col_fun = boxcol, 
		       title = boxtitle),
                Legend(col_fun = outlcol, 
		       title = outltitle)
        )

        hm <- hm.cell(boxdat,
                right_annotation = ha,
                cell.w = .15,
                cell.h = .15,
                show_heatmap_legend = F,
                col = boxcol
        )

        dir.pdf(out, path, 
		height = 24, 
		width = 10 + ncol(boxdat), 
		append.date = F)
        draw(hm, annotation_legend_list = lgd)
        dev.off()
}

clustparam <- function(dat, clusts, out,
		       logfc.cutoff = 0.5, 
		       fdr.cutoff = 0.05){
	clustdat <- split(dat, clusts)
	test <- test.clust.params(dat,clustdat)

	# select fields from output
	mudat <- sapply(test, '[', "mu", )
	fcdat <- sapply(test, '[', "FC", )
	p.t <- sapply(test, '[', 't', )
	p.u <- sapply(test, '[', 'u', )

	# convert p-values to FDR values
	fdr.t <- apply(p.t, 2, 
		       function(x) p.adjust(unlist(x)))
	fdr.u <- apply(p.u, 2, 
		       function(x) p.adjust(unlist(x)))

	#fdr.t <- p.adjust(p.t)
	#fdr.u <- p.adjust(p.u)

	log.fc <- log2(fcdat)
	log.fc[!is.finite(log.fc)] <- 0

	log.t <- -log10(fdr.t)
	log.u <- -log10(fdr.u)

	sel <- abs(log.fc) > logfc.cutoff
	sel.t <- fdr.t <= fdr.cutoff & sel
	sel.u <- fdr.u <= fdr.cutoff & sel

	rsel.t <- apply(sel.t,1,any)
	rsel.u <- apply(sel.u,1,any)

	# define color scale
	clustdat <- lapply(clustdat, 
			   function(x) t(x)[row.names(log.fc),])
	clustdat.t <- lapply(clustdat,`[`,rsel.t,) 
	clustdat.u <- lapply(clustdat,`[`,rsel.u,) 
	# col.z creates a color scale centered on 0 ranging 
	# from the 0.01 to 0.99 quantiles of the input data
	if(length(unique(clusts)) > 1){
		box.heatmap(m[rsel.t,], 
			    clustdat.t,
			    log.fc[rsel.t,], 
			    log.t[rsel.t,], 
			    "t.boxplot", out) 

		box.heatmap(m[rsel.u,], 
			    clustdat.u,
			    log.fc[rsel.u,], 
			    log.u[rsel.u,], 
			    "u.boxplot", out) 

		writepdf({
			dotplot(log.fc[rsel.t,], 
				log.t[rsel.t,],
				mat.name = "log2(FC)",
				row_title_rot = 0)
		}, "t.fc.pdf", out, height = 24)

		writepdf({
			dotplot(log.fc[rsel.u,], 
				log.u[rsel.u,],
				mat.name = "log2(FC)",
				row_title_rot = 0)
		}, "u.fc.pdf", out, height = 24)
	
		writepdf({
			dotplot(mudat[rsel.t,], 
				log.t[rsel.t,],
				mat.name = "mean",
				row_title_rot = 0)
		}, "t.pdf", out, height = 24)

		writepdf({
			dotplot(mudat[rsel.u,], 
				log.u[rsel.u,],
				mat.name = "mean",
				row_title_rot = 0)
		}, "u.pdf", out, height = 24)
	}
}


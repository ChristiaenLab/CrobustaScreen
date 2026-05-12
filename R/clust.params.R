source("R/dotplot.R")

clust.t.test <- function(x, clust, dat) {
        mu <- mean(dat[, x], na.rm = T)
        mu_cl <- mean(clust[, x], na.rm = T)
        FC <- mu_cl / mu
        if(FC != 1 & length(unique(clust[, x])) > 1) {
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

test.clust.params <- function(dat, clustdat) {
# test for significance of each feature in each cluster
    clusttest <- function(clust) {
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
          outltitle = "-log10(FDR)", ylim = c(-1, 1), 
		  omit = NULL, ...) {
		if(!is.null(omit)) m <- m[, !omit]
        # clip clustdat to ylim
        clustdat <- lapply(clustdat, function(x) {
                x[x < ylim[1]] <- ylim[1]
                x[x > ylim[2]] <- ylim[2]
                return(x)
        })

    # This scale will be used for the fold changes
    # between average feature values in each cluster
    # from the background.
        boxcol <- col.z(boxdat, .05, 0)

        # creates a color scale from 0 to 2
        # This will give a log10 scale for FDR values
    # between 1 and 0.01
        outlcol <- colorRamp2(c(0, 2), c("white", "black"))

        # show feature values within clusters as boxplots
    # color-coded by FDR value
        getAnn <- function(x) {
                m <- clustdat[[x]]
                fc.cols <- boxcol(boxdat[, x])
                fdr.cols <- outlcol(outldat[, x])
                return(anno_boxplot(
                        m,
                        which = "row",
                        width = unit(1, "in"),
                        box_width = 0.9,
                        ylim=ylim,
                        gp = gpar(
                                fill = fc.cols,
                                col = fdr.cols
                        ), ...
                ))
        }

        # apply annotation function to each cluster
        ha <- lapply(colnames(boxdat), getAnn)
        #names(ha) <- paste("cluster", colnames(boxdat))
        names(ha) <- colnames(boxdat)

        # bind annotation into single object
        ha <- do.call(HeatmapAnnotation,
              append(ha, list(which = "row")))

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

clustparam <- function(m, clusts, path,
               logfc.cutoff = 0.5,
               fdr.cutoff = 0.05,
               subset = NULL,
               omit = NULL, ...) {
    clustdat <- split(m, clusts)
    test <- test.clust.params(m, clustdat)

    if(!is.null(subset)) {
        test <- test[subset]
    }

    # select fields from output
    mudat <- sapply(test, "[", "mu", )
    fcdat <- sapply(test, "[", "FC", )
    p.t <- sapply(test, "[", "t", )
    p.u <- sapply(test, "[", "u", )

    # convert p-values to FDR values
    fdr.t <- apply(p.t, 2,
               function(x) p.adjust(unlist(x)))
    fdr.u <- apply(p.u, 2,
               function(x) p.adjust(unlist(x)))

    dir.csv(mudat, "mean", path, append.date = F)
    dir.csv(fcdat, "FC", path, append.date = F)
    dir.csv(fdr.t, "FDR_t", path, append.date = F)
    dir.csv(fdr.u, "FDR_u", path, append.date = F)

    log.fc <- log2(fcdat)
    log.fc[!is.finite(log.fc)] <- 0

    log.t <- -log10(fdr.t)
    log.u <- -log10(fdr.u)

    sel <- abs(log.fc) > logfc.cutoff
    sel.t <- fdr.t <= fdr.cutoff & sel
    sel.u <- fdr.u <= fdr.cutoff & sel

    rsel.t <- apply(sel.t, 1, any)
    rsel.u <- apply(sel.u, 1, any)

    # define color scale
    clustdat <- lapply(clustdat,
               function(x) t(x)[row.names(log.fc),])
    clustdat.t <- lapply(clustdat,`[`, rsel.t,)
    clustdat.u <- lapply(clustdat,`[`, rsel.u,)

    if(length(unique(clusts)) > 1) {
        box.heatmap(m[rsel.t,],
                clustdat.t,
                log.fc[rsel.t,],
                log.t[rsel.t,],
                "t.boxplot", path,
				omit = omit, ...)

        box.heatmap(m[rsel.u,],
                clustdat.u,
                log.fc[rsel.u,],
                log.u[rsel.u,],
                "u.boxplot", path,
				omit = omit, ...)

        dotplot(log.fc[rsel.t,],
                log.t[rsel.t,],
                mat.name = "log2(FC)",
                row_title_rot = 0,
                filename = "t.fc.pdf", path = path)

        dotplot(log.fc[rsel.u,],
                log.u[rsel.u,],
                mat.name = "log2(FC)",
                row_title_rot = 0,
                filename = "u.fc.pdf", path = path)
   
        dotplot(mudat[rsel.t,],
                log.t[rsel.t,],
                mat.name = "mean",
                row_title_rot = 0,
                filename = "t.pdf", path = path)

        dotplot(mudat[rsel.u,],
                log.u[rsel.u,],
                mat.name = "mean",
                row_title_rot = 0,
                filename = "u.pdf", path = path)
    }
}

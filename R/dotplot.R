source("R/dirfns.R")
source("R/heatmapfns.R")

writepdf <- function(expr, file, out='.', append.date=F,...) {
    out <- mkdate(file, ext='', path=out, append.date=append.date)
    pdf(out,...)
    tryCatch(expr, finally=dev.off())
}

#' accepts the results of an enrichment test applied to each cell in a matrix
#' and writes a dotplot of the results

#' @param mat A matrix of values shown by the dot color.
#' @param outl A matrix of values shown by the dot outlline.
#' @param size A matrix of values shown by the dot size.
#' @param col.mat A color scale for \code{mat}.
#' @param col.outl A color scale for \code{outlline}.
#' @param scale A vector of length 2 giving the min and max values to scale the size of the dots.
#' @param cell.dim The width & height of each heatmap cell in inches.
#' @param ... Additional arguments to \code{hm.cell()}.
#' @export
#' @importFrom grid gpar unit grid.points
hmdot.outl <- function(
    mat, outl, size, 
    col.mat, col.outl, scale, size.breaks,
    mat.name="log2(OR)", 
    outl.name="size", 
    size.name="-log10(FDR)", 
        cell.dim=.15,
    filename = NULL, path = '.', append.date = F, # Added for dir.hm
    ...
) {
    require(ComplexHeatmap)
    mat[is.na(mat)] <- 0
    mat[mat==Inf] <- max(mat[is.finite(mat)])
    mat[mat==-Inf] <- min(mat[is.finite(mat)])

    cexfn <- function(x) unit((1.2 * x / max(size)) * 
                  cell.dim, 'in')
    cellfn <- function(j, i, x, y, width, height, fill) {
            grid.points(
        x = x, y = y, 
        size = cexfn(size[i, j]),
        pch = 16,
                gp = gpar(
            col = col.mat(mat[i, j]) 
        )
        )
            grid.points(
        x = x, y = y, 
        size=cexfn(size[i, j]),
        pch=1,
                gp = gpar(
            col = col.outl(outl[i, j])
        )
        )
        }

    lgd <- list(
        Legend(col_fun = col.mat, title = mat.name),
        Legend(col_fun = col.outl, title = outl.name),
        Legend( title=size.name,
            at=size.breaks,
            type='points',
            background=0,
            pch=16,
            size=unit(sapply(size.breaks, cexfn),'in'),
            legend_gp=gpar(col=1, fill=0)
        )
    )

    hm <- hm.cell(
        mat,
        cell_fun=cellfn,
        rect_gp = gpar(type = "none"),
        cell.w=cell.dim,
        cell.h=cell.dim,
        show_heatmap_legend=F,
        ...
    )

    # Use dir.hm to save with automatic sizing
    if (!is.null(filename)) {
        dir.hm(hm, filename, path = path, annotation_legend_list = lgd, append.date = append.date)
    }
    return(invisible(hm)) # Return hm invisibly for potential further manipulation
}

hmdot <- function(
    mat, size, 
    col.mat, scale, size.breaks,
    mat.name="log2(OR)", 
    size.name="-log10(FDR)", 
        cell.dim=.15,
    filename = NULL, path = '.', append.date = F, # Added for dir.hm
    ...
) {
    require(ComplexHeatmap)
    mat[is.na(mat)] <- 0
    mat[mat==Inf] <- max(mat[is.finite(mat)])
    mat[mat==-Inf] <- min(mat[is.finite(mat)])

    cexfn <- function(x) unit((1.2 * x / max(size)) * 
                  cell.dim, 'in')
    cellfn <- function(j, i, x, y, width, height, fill) {
            grid.points(
        x = x, y = y, 
        size = cexfn(size[i, j]),
        pch = 16,
                gp = gpar(
            col = col.mat(mat[i, j]) 
        )
        )
        }

    lgd <- list(
        Legend(col_fun = col.mat, title = mat.name),
        Legend( title=size.name,
            at=size.breaks,
            type='points',
            background=0,
            pch=16,
            size=unit(sapply(size.breaks, cexfn),'in'),
            legend_gp=gpar(col=1, fill=0)
        )
    )

    hm <- hm.cell(
        mat,
        cell_fun=cellfn,
        rect_gp = gpar(type = "none"),
        cell.w=cell.dim,
        cell.h=cell.dim,
        show_heatmap_legend=F,
        ...
    )

    # Use dir.hm to save with automatic sizing
    if (!is.null(filename)) {
        dir.hm(hm, filename, path = path, annotation_legend_list = lgd, append.date = append.date)
    }
    return(invisible(hm)) # Return hm invisibly
}


#' @param dot dot color matrix
#' @param size dot size matrix
#' @param outline dot outline color matrix
#' @param sizelim Maximum value on the size scale. Values above \code{sizelim} are set to \code{sizelim}.
#' @param ... Additional arguments to \code{hmdot()}.
#' @export
dotplot.outl <- function(dot, size, outline, 
            sizelim = -log10(0.01), 
            outl.name = 'size', 
            size.name = '-log10(FDR)', ...) {
    size[which(size > sizelim)] <- sizelim
    #logFDR <- -log10(size)
    dot[is.na(dot)] <- 0
    dot[dot == Inf] <- sizelim
    dot[dot == -Inf] <- 0
    size[!is.finite(size)] <- 0

    col.dot <- col.z(dot)

    outline[!is.finite(outline)] <- 0
    outlinescale <- c(0, quantile(as.matrix(outline)[
                       as.matrix(outline) != 0 
                      ], 0.95))
    col.outl <- colorRamp2(outlinescale, c('white', 'black'))

    size.breaks <- seq(0, sizelim, length.out = 6)

    hmdot.outl(dot, 
           outline, 
           size, 
           col.mat = col.dot, 
           col.outl = col.outl, 
           scale = outlinescale, 
           size.breaks = size.breaks,
           outl.name = outl.name,
           size.name = size.name, ...)
}

dotplot <- function(dot, size,
            sizelim = -log10(0.01), 
            size.name = '-log10(FDR)', ...) {
    size[which(size > sizelim)] <- sizelim
    dot[is.na(dot)] <- 0
    dot[dot == Inf] <- sizelim
    dot[dot == -Inf] <- 0
    size[!is.finite(size)] <- 0

    col.dot <- col.z(dot)

    size.breaks <- seq(0, sizelim, length.out = 6)

    hmdot(dot, 
          size, 
          col.mat = col.dot, 
          size.breaks = size.breaks,
          size.name = size.name, ...)
}

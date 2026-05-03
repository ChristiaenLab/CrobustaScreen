library(circlize)

#' Color scale for a specified quantile. This scale is intended for heatmaps containing negative and positive values, so the range is set to \code{c(quant,1-quant)}.
#'
#' @param x A numeric matrix or vector.
#' @param quant The quantile to be used as the scale limits.
#' @param mid The midpoint of the scale.
#' @return A \code{colorRamp2} scale which can be passed to \code{Heatmap()}.
#' @importFrom circlize colorRamp2
#' @export
col.z <- function(x,quant=.01, mid=0, cols=c('blue', 'white', 'red')) {
	breaks <- c(quantile(x, quant, na.rm=T),
		    mid,
		    quantile(x, 1-quant, na.rm=T))
	colorRamp2(breaks, cols)
}


#' Color scale for a specified quantile. This scale is intended for heatmaps containing only positive values, so the range is set to \code{c(0,1-quant)}.
#'
#' @param x A numeric matrix or vector.
#' @param quant The quantile to be used as the upper limit.
#' @param cols The colors used for the color scale.
#' @return A \code{colorRamp2} scale which can be passed to \code{Heatmap()}.
#' @importFrom circlize colorRamp2
#' @export
col.abs <- function(x, quant=.05, min=0, cols=c('white', 'black')){
	breaks <- c(min, quantile(x[x!=0], 1 - quant, na.rm=T))
	breaks <- seq(breaks[1], breaks[2], length.out=length(cols))
	colorRamp2(breaks, cols)
}

#' Color scale for categorical data.
#'
#' @param cond A vector that can be coerced to a factor.
#' @param colfn A function that returns a color map for each level in \code{cond}.
#' @param ... Additional arguments to \code{colfn}.
#' @return A named vector of colors corresponding to the levels of \code{cond}.
#' @export
cond.col <- function(cond, colfn=rainbow,...){
	cond <- as.factor(cond)
	cols <- colfn(length(levels(cond)), ...)
	names(cols) <- levels(cond)
	return(cols)
}

#' Creates a color scale for the levels in a vector, then returns a vector assigning a color to each element of the input vector.
#'
#' @param cond A vector that can be coerced to a factor.
#' @param colfn A function that returns a color map for each level in \code{cond}.
#' @param ... Additional arguments to \code{colfn}.
#' @return A vector of colors corresponding to the elements of \code{cond}.
#' @export
cond.col.vec <- function(cond, colfn=rainbow,...) {
	cond.col(cond, colfn=colfn, ...)[as.numeric(as.factor(cond))]
}

writepdf <- function(expr, file, out='.',...){
    pdf(paste0(out, '/', file), ...)
    tryCatch(expr, finally=dev.off())
}
        
#' Wrapper for \code{Heatmap()} which allows specifying cell dimensions and resizing the heatmap accordingly.
#'
#' @param x A numeric matrix to be plotted.
#' @param ... Additional arguments to \code{Heatmap()}.
#' @param cell.h The cell height.
#' @param cell.w The cell width.
#' @param height The heatmap height. Ignored if \code{cell.h} is specified.
#' @param width The heatmap width. Ignored if \code{cell.w} is specified.
#' @param units The unit scale to be used for \code{cell.h} and \code{cell.w}.
#' @return A ComplexHeatmap.
#' @import ComplexHeatmap
#' @export
hm.cell <- function(
		x,...,
		cell.h=NULL, cell.w=NULL,
		height=NULL, width=NULL,
		# heatmap_height=NULL, 
		# heatmap_width=NULL,
		units='in'
){
	if(!is.null(cell.h)) height <- unit(nrow(x) * cell.h, units)
	if(!is.null(cell.w)) width <- unit(ncol(x) * cell.w, units)
	return(Heatmap(x, ..., height=height, width=width))
}

hm.quant <- function(x, ..., quant=.01, mid=0, cols=c('blue', 'white', 'red')) {
	col <- col.z(x, quant, mid, cols)
	hm.cell(x, ...,col=col)
}

hm.abs <- function(x, ..., quant=.05, min=0, cols=c('white', 'black')) {
	col <- col.abs(x, quant, min, cols)
	hm.cell(x, ...,col=col)
}

#' Remove columns with zero variance
#'
#' @param x A data.frame or matrix.
#' @return The input with zero-variance columns removed.
rm.zero.var <- function(x) {
	vars <- apply(x, 2, var, na.rm=TRUE)
	x[, vars > 0, drop=FALSE]
}

#' Correlation heatmap
#'
#' @param x A data.frame or matrix. Columns with zero variance are removed.
#' @param method The correlation method passed to \code{cor()}.
#' @param ... Additional arguments to \code{hm.cell()}.
#' @param cols The colors for the correlation scale (negative, zero, positive).
#' @return A ComplexHeatmap showing correlations between all columns of \code{x}.
#' @export
cor.hm <- function(x, method='pearson', ..., cols=c('blue', 'white', 'red')) {
	x <- rm.zero.var(x)
	C <- cor(x, method=method, use='pairwise.complete.obs')
	col <- colorRamp2(c(-1, 0, 1), cols)
	hm.cell(C, ..., col=col)
}

#' Cross-correlation heatmap between two tables
#'
#' @param x A data.frame or matrix. Columns with zero variance are removed.
#' @param y A data.frame or matrix with the same number of rows as \code{x}. Columns with zero variance are removed.
#' @param method The correlation method passed to \code{cor()}.
#' @param ... Additional arguments to \code{hm.cell()}.
#' @param cols The colors for the correlation scale (negative, zero, positive).
#' @return A ComplexHeatmap showing correlations between columns of \code{x} and columns of \code{y}.
#' @export
xycor.hm <- function(x, y, method='pearson', ..., cols=c('blue', 'white', 'red')) {
	x <- rm.zero.var(x)
	y <- rm.zero.var(y)
	C <- cor(x, y, method=method, use='pairwise.complete.obs')
	col <- colorRamp2(c(-1, 0, 1), cols)
	hm.cell(C, ..., col=col)
}

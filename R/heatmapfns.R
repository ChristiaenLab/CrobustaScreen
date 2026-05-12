library(circlize)
library(ComplexHeatmap)

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
		cell.h=NULL,cell.w=NULL,
		height=NULL,width=NULL,
		# heatmap_height=NULL, 
		# heatmap_width=NULL,
		units='in'
){
	if(!is.null(cell.h)) height <- unit(nrow(x)*cell.h,units)
	if(!is.null(cell.w)) width <- unit(ncol(x)*cell.w,units)
	return(Heatmap(x,...,height=height,width=width))
}

#' Heatmap with color scale based on quantiles.
#' @param x A numeric matrix.
#' @param filename Output filename.
#' @param path Output path.
#' @param ... Additional arguments to Heatmap.
#' @export
quantHeatmap <- function(x, filename = NULL, path = '.', 
                         name = "value",
                         show_row_names = TRUE, show_column_names = TRUE,
                         show_row_dend = TRUE, show_column_dend = TRUE,
                         ..., 
                         append.date = F){
	require(ComplexHeatmap)
	col <- col.z(x)
	hm <- hm.cell(x, col = col, 
                  name = name,
                  show_row_names = show_row_names, show_column_names = show_column_names,
                  show_row_dend = show_row_dend, show_column_dend = show_column_dend,
                  ...)
	
	if(!is.null(filename)){
		dir.hm(hm, filename, path = path, append.date = append.date)
	}
	return(hm)
}

#' Color scale for a specified quantile. This scale is intended for heatmaps containing negative and positive values, so the range is set to \code{c(quant,1-quant)}.
#'
#' @param x A numeric matrix or vector.
#' @param quant The quantile to be used as the scale limits.
#' @param mid The midpoint of the scale.
#' @return A \code{colorRamp2} scale which can be passed to \code{Heatmap()}.
#' @importFrom circlize colorRamp2
#' @export
col.z <- function(x,quant=.01, mid=0, cols=c('blue', 'white', 'red')) {
	lower <- quantile(x, quant, na.rm = T)
	upper <- quantile(x, 1 - quant, na.rm = T)
	if(upper > 0) {
		if(lower < 0) {
			breaks <- c(lower, mid, upper)
			colorRamp2(breaks, cols)
		} else {
			breaks <- c(lower, upper)
			colorRamp2(breaks, cols[-1])
		}
	} else {
		breaks <- c(lower, upper)
		colorRamp2(breaks, cols[-2])
	}
}


#' Color scale for a specified quantile. This scale is intended for heatmaps containing only positive values, so the range is set to \code{c(0,1-quant)}.
#'
#' @param x A numeric matrix or vector.
#' @param quant The quantile to be used as the upper limit.
#' @param cols The colors used for the color scale.
#' @return A \code{colorRamp2} scale which can be passed to \code{Heatmap()}.
#' @importFrom circlize colorRamp2
#' @export
col.abs <- function(x,quant=.05, cols=c('white','black')){
	breaks <- c(0,quantile(x[x!=0],1-quant, na.rm=T))
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
cond.col <- function(cond,colfn=rainbow,...){
	cond <- as.factor(cond)
	cols <- colfn(length(levels(cond)),...)
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
cond.col.vec <- function(cond,colfn=rainbow,...) {
	cond.col(cond,colfn=colfn,...)[as.numeric(as.factor(cond))]
}


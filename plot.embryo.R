source("imageFns.R")
require(purrr)

get.axes <- function(x,prefix="Nucleus"){
	sapply(c("A","B","C"),
	       compose(partial(coord,x, prefix), 
		       partial(paste0,"Ellipsoid.Axis.")), 
	       simplify=F)
}

cell.axes <- function(i,ax) do.call(rbind, 
				  lapply(ax,'[',1,,drop=T))

flatten.ellipsoid <- function(ax){
	dists <- apply(ax[,-3],1,
		       compose(dist,
			       partial(rbind,c(0,0))))
	return(ax[order(abs(dists),
			decreasing=T)[-3],-3])
}

axes <- function(x,prefix="Nucleus",i){
	compose(flatten.ellipsoid,
		partial(cell.axes,i),
		get.axes)(x,prefix)
}

ellipses <- function(x,prefix="Nucleus"){
	sapply(1:nrow(x),partial(axes,x,prefix),simplify=F)
}

lim <- function(pos,ax,i){
	ax <- apply(do.call(rbind,ax),2,unlist)
	return(c(min(pos[,i])-min(ax[,i]),
	  max(pos[,i])+max(ax[,i])))
}

plot.embryo <- function(cells,surfaces=NULL){
	nucl.pos <- ncoord(cells)[,-3]
	cell.pos <- ccoord(cells)[,-3]
	surface.pos <- coord(surfaces,"")[,-3]
	nucl.ax <- ellipses(cells)
	cell.ax <- ellipses(cells,"Cell")
	surface.ax <- ellipses(surfaces,'')

	plot(NULL,
	     xlim=lim(surface.pos,surface.ax,1),
	     ylim=lim(surface.pos,surface.ax,2))

	fn <- function(pos,ax) {
		function(...){
			sapply(1:nrow(pos),
			       function(i) ellipse(pos[i,],
						   ax[[i]],
						   ...))
		}
	}

	fn(nucl.pos,nucl.ax)(col='red')
	fn(cell.pos,cell.ax)(col='forestgreen')
	fn(surface.pos,surface.ax)()
}

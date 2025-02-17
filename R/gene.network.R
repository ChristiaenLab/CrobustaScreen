source("R/dirfns.R")
source("R/pois.R")

dotPois <- function(pois, out, append.date = T){
	require(moreComplexHeatmap)
	#require(dirfns)

	dotPscale(
		pois$log2OR, 
		pois$p, 
		pois$count, 
		file = 'conditionEdgePois', 
		path = out, 
		row_title_rot = 0,
		append.date = append.date
	)

	dir.csv(pois$log2OR, 'conditionPoisLog2OR', 
		out, append.date = append.date)
	dir.csv(pois$p, 'conditionPoisFDR', 
		out, append.date = append.date)
	dir.csv(pois$count, 'conditionPoisCount', 
		out, append.date = append.date)
}

poisGraph <- function(pois, cutoff = 0.05, up = T, mode = 'max'){
	require(igraph)
	require(circlize)

	edgeval <- pois$log2OR
	edgeval[which(!is.finite(edgeval))] <- 0
	edgeval[which(pois$p > cutoff)] <- 0
	if(up){
		edgeval[which(edgeval < 0)] <- 0
	}
	g <- graph_from_adjacency_matrix(edgeval,
					 mode, T, F)
	if(up) { 
		colfn <- colorRamp2(c(0, max(E(g)$weight)),
				   c('white', 'red')) 
	}else colfn <- colorRamp2(c(min(E(g)$weight), 0,
				    max(E(g)$weight)),
				  c('blue', 'white', 'red'))
	E(g)$edge.color <- colfn(E(g)$weight)

	return(g)
}

networkPois <- function(g, file, out = '.',
			colfn = col.z(E(g)$weight),
			fn = plotNetwork,
			title = 'log2OR',
			layout = layout_nicely, ...){
	require(igraph)
	#require(dirfns)

	E(g)$edge.color <- colfn(E(g)$weight)

	lgd <- seq(round(quantile(E(g)$weight, 0.99)),
		   min(0, round(quantile(E(g)$weight,
					0.01))),
		   length.out = 6)
	fn(g, file, out, colfn, lgd, title, layout, ...)
}

plotNetwork <- function(g, file, out = '.', colfn, lgd,
			title = 'log2OR', layout = layout_nicely,
			vertex.shape = 'none',
			vertex.size = 10, ...,
			append.date = T){
	tmp <- do.call(rbind, strsplit(as_ids(E(g)), '\\|'))
	curved <- sapply(1:nrow(tmp), function(x) { 
		y <- which(tmp[, 1] == tmp[x, 2] & tmp[, 2] == tmp[x, 1]) 
		if(length(y) > 0) return(c(x, y))
	})

	curves <- curve_multiple(g)
	curves[unique(unlist(curved))] <- 0.3
	#         lgd <- seq(lgd[1], 0, length.out = 6)

	#         ldist <- distances(g, weights = abs(E(g)$weight))
	#         ldist[is.finite(ldist) & ldist > 0] <- 1
	#         l <- layout_with_mds(g, dist = ldist, dim = 2, 
	#                         options = arpack_defaults)

	dir.pdf(file, out, append.date = append.date)
	E(g)$weight <- 1
	plot(g, vertex.shape = vertex.shape,
	     vertex.size = vertex.size,
	     layout = layout,
	     edge.curved = curves,
	     #autocurve = T,
	     edge.arrow.size = 0.5,
	     #layout = l,
	     edge.color = E(g)$edge.color, #colfn(E(g)$weight),
	     ...)
	legend('topleft', as.character(lgd), fill = colfn(lgd), title = title)
	dev.off()

}

plotNetworkCircle <- function(g, file, out = '.', colfn, lgd, 
			      title = 'log2OR', 
			      layout = layout.circle, 
			      append.date = T){
	#require(dirfns)

	dir.pdf(file, out, append.date = append.date)
	plot(g, vertex.shape = 'none', 
	     layout = layout.circle, 
	     edge.color = E(g)$edge.color)
	legend('topleft', as.character(lgd), 
	       fill = colfn(lgd), title = title)
	dev.off()
}

enrichCond <- function(cond, dists, out, layout = layout_nicely){
	require(moreComplexHeatmap)
	require(igraph)

	#number of embryos per condition
	ncond <- table(cond)

	pois <- getPois(cond, dists)
	dotPois(pois, out)

	g <- poisGraph(pois, up = F)
	colfn <- col.z(E(g)$weight)
	networkPois(g, 'conditionEdgeNetwork', out, colfn, layout = layout.circle, vertex.shape = 'none')

	g.up <- poisGraph(pois)
	networkPois(g.up, 'conditionEdgeNetworkUp',
		    out, colfn, layout = layout)
}



wak <- function(g){
	require(igraph)

	G <- as.matrix(as_adjacency_matrix(g))
	G / apply(G,2,sum)
}

diffuse <- function(g, X){
	X <- as.matrix(X)
	G <- wak(g)
	G %*% X
}

mse <- function(X, Xhat){
	mean((X - Xhat) ^ 2)
}

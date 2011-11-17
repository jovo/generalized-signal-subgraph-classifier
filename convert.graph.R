convert.graph <- function(graph){
	## Data cannot have NA's
	if (!all(complete.cases(c(graph)))) stop("Data cannot have NA's")
	## get class of data so can handle multiple inputs
	cg <- class(graph)
	if (!cg %in% c("array", "data.frame", "matrix")) stop("Data must be array/data.frame/matrix")
	## see if data is stored as array
	if (cg == "array"){
		dag <- dim(graph)
		## 3rd dimension is always person/unit
		if (length(dag) > 3) stop("Too Many dimensions, max = 3")
		## If taking Upper.triangular
		if (upper) {
			## need square matrix
			if (dag[1] != dag[2]) stop("Upper Selected but not square matrix")
			graph <- t(apply(graph, 3, function(x) return(x[upper.tri(x)])))
		} else {
			## converting m x n x p array to p x mn matrix
			graph <- t(matrix(c(graph), nrow=dag[1]* dag[2], ncol = dag[3]))
		}
	}
	## convert data.frame (one column is labels)
	if (cg == "data.frame") {
		## if ylabels are a column in data.frame, then pull it out
		ng <- colnames(graph)
		classes <- unique(sapply(ng, function(x) class(graph[,x])))
		if (!all(classes %in% c("numeric", "logical", "integer"))) stop("Some Column has wrong class")
		graph <- as.matrix(graph[, ng[!ng %in% "y"]])
	}
	## Last option is data is matrix
	return(graph)
	
}



getpvals <- function(x, y, workspace, weight, ...){
	if (weight){
		## Kruskal - Wallis test for pvals for weighted
		return(kruskal.test(x=x, g=y, ...)$p.value)
	} else {
		tab <- table(x, ylabels)
		expected <- rowSums(tab) %*% t(colSums(tab))/sum(tab)
		if (all(expected > 5) && (nrow(tab) > 2 || ncol(tab) > 2)) {
			return(chisq.test(x=x, y=y)$p.value)
			### stopped here
		} else return(fisher.test(x=x, y=y, alternative=alternative, conf.int = FALSE, workspace=workspace)$p.value)
	}	

	
}


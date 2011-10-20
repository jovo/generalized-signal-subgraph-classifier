rm(list=ls())
source("convert.graph.R")
#### Test data
graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 1, 0.5)
#graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
#graph <- data.frame(cbind(graph, ylabels=ylabels))
upper <- FALSE
alternative = "two.sided"
keep.pct=0.1

#binary signal subgraph
bssgraph <- function(graph, ylabels=NULL, upper=FALSE, alternative = "two.sided", keep=0.5, workspace=400000, ...){
	## taken from fisher.test <- check if alternative is in acceptable range
	alternative <- char.expand(alternative, c("two.sided", "less", "greater"))
	if (length(alternative) > 1L || is.na(alternative)) stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
	
	dag <- dim(graph)
	cg <- class(graph)
	## see if data is stored as array
	## if ylabels are a column in data.frame, then pull it out
	if (cg == "data.frame" & ("ylabels" %in% colnames(graph)))
		ylabels <- graph$ylabels	
	if (cg == "matrix" && ("ylabels" %in% colnames(graph))) {
		ylabels <- graph[,"ylabels"]
		which.y <- which(colnames(graph) == "ylabels")
		graph <- graph[,-which.y]
	}
	
	# convert into correct format
	graph <- convert.graph(graph)

	## get levels of adjacencies
	adj <- sort(unique(c(graph)))
	nadj <- length(unique(c(graph)))
	
	if (nadj > 2) stop("Too many labels - graph is not in {0, 1}, maybe run wssgraph?")
	
	## Get number of groups
	ngroups <- length(unique(ylabels))
	if (alternative != "two.sided" & ngroups > 2) {
		print("For groups/labels > 2, alternative is two sided")
		alternative <- "two.sided"
	}

	levs <- sort(unique(ylabels))
	if (ngroups == 2) {
		if (class(levs) == "character") ylabels <- as.numeric(factor(ylabels, levels=levs))-1
		if (class(levs) == "factor") ylabels <- as.numeric(ylabels, levels=levs)-1
		if (!all(as.numeric(levs)==c(0, 1))) stop("Label problems")
	}
	dg <- dim(graph)
	## get number of subjects/vertices
	nvert <- dg[2]
	nsubj <- dg[1]
	## rows are units/subjects, columns are vertices
	if (length(dg) > 2) stop("Problem with dimensions")
	
	## get p.value for the graph vertices
	
	pvals <- apply(graph, 2, function(x) getpvals(x, ylabels=ylabels, workspace=workspace))
	## keep percentage of the graph (need to CV it probably)
	nkeep <- floor(nvert*keep.pct)
	
	## ssgraph: indices of the signal subgraph
	ssgraph <- order(pvals)[1:nkeep]
	
	if (cg == "array") {
		sig.graph <- matrix(0, nrow=dag[1], ncol=dag[2])
		## make sure to put the ssgraph back in the right spot if upper is on
		if (upper){
			up.tri <- upper.tri(sig.graph)
			sig.graph[up.tri][ssgraph] <- 1
		} else sig.graph[ssgraph] <- 1
	} else {
		sig.graph <- rep(0, nvert)
		sig.graph[ssgraph] <- 1
	}
	
	## only need product over signal subgraph
	## tgraph is total graph
	tgraph <- graph
	graph <- graph[, ssgraph]


	## Proportion matrices for each group
	## Each column is a different group
	p.mat <- matrix(nrow=nkeep, ncol=ngroups)
	priors <- rep(NA, length=ngroups)	
	for (ig in 1:ngroups){
		## grab the level of ylabels
		ilev <- levs[ig]
		## get indicator if each person is in that group
		ingroup <- ylabels == ilev
		# priors = P(Y = y)
		priors[ig] <- mean(ingroup)
		## Subset the graphs to only that group
		sset <- graph[ingroup,]
		## get the mean of the adjacency
		## p.mat <- p_{u,v | Y = y}
		p.mat[, ig] <- colMeans(sset)
	}

#getprob <- function(graph, priors) {
	# probs 
	probs <- matrix(nrow=nsubj, ncol=ngroups)
	colnames(probs) <- levs
	
	for (ig in 1:ngroups){
#		for (ikeep in 1:nkeep){
			ps <- dbinom(graph, 1, prob=p.mat[,ig])
			probs[, ig] <- apply(ps, 1, prod)*priors[ig]
#		}
	}
	## divide by total to get actual probability and not just numerator P(G = g | Y = y)
	## apply(probs, 1, sum) = P(G=g|Y=1)P(Y=1)+ ... + P(G=g|Y=y_n)P(Y=y_n)
	probs <- probs / apply(probs, 1, sum)
	
	preds <- apply(probs, 1, function(x) which(x == max(x)))
	if (class(preds) != "integer") stop("Problem with Prediction")
	preds <- levs[preds]
	print(table(preds, ylabels))
	print(mean(preds == ylabels))
	
	return(list(ssgraph=sig.graph, priors=priors, prob=probs, pred=preds, ylabels=ylabels))
#}

}

graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 1, 0.5)
x <- bssgraph(agraph, ylabels=ylabels)
x <- bssgraph(agraph, ylabels=ylabels, alternative="less")
x <- bssgraph(agraph, ylabels=ylabels, alternative="great")
graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
x <- bssgraph(graph, ylabels=ylabels)

x <- bssgraph(graph, ylabels=ylabels, alternative="less")
x <- bssgraph(graph, ylabels=ylabels, alternative="great")

graph <- data.frame(cbind(graph, ylabels=ylabels))

x <- bssgraph(graph, ylabels=ylabels)
x <- bssgraph(graph, ylabels=ylabels, alternative="less")
x <- bssgraph(graph, ylabels=ylabels, alternative="great")

graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 5, 0.5)
x <- bssgraph(agraph, ylabels=ylabels, workspace=10000000)
x <- bssgraph(agraph, ylabels=ylabels, alternative="less")
x <- bssgraph(agraph, ylabels=ylabels, alternative="great")


setwd("~/Dropbox/CTR/Stewart/Signal_Subgraph/")
rm(list=ls())
source("convert.graph.R")
#### Test data
Cl <- cor(longley)
graph <- agraph <- array(runif(1000 * 100, min=-1, max=1), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 5, 0.5)
alternative = "two.sided"
upper = TRUE
keep.pct=0.1
corr=TRUE


#binary signal subgraph
wssgraph <- function(graph, ylabels=NULL, upper=FALSE, keep=NULL, corr=FALSE, ...){
	
	
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
	
	# convert into correct format - rows are units/subjects, columns are nodes
	graph <- convert.graph(graph)

	adj <- sort(unique(c(graph)))
	nadj <- length(unique(c(graph)))
	
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
	pvals <- apply(graph, 2, function(x) getwpvals(x, ylabels=ylabels))
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

	if (corr) graph <- (graph+1)/2

	## Proportion matrices for each group
	## Each column is a different group
	v.mat <- sd.mat <- m.mat <- matrix(nrow=nkeep, ncol=ngroups)
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
		m.mat[, ig] <- apply(sset, 2, mean)
		sd.mat[, ig] <- apply(sset, 2, sd)
		v.mat[, ig] <- apply(sset, 2, var)
	}

#getprob <- function(graph, priors) {
	# probs 
	probs <- matrix(nrow=nsubj, ncol=ngroups)
	colnames(probs) <- levs
	
	if (corr==FALSE) {
		for (ig in 1:ngroups){
			ps <- matrix(nrow=nsubj, ncol=nkeep)
			for (ikeep in 1:nkeep){
				#ps is the probability given the mean/sd of each group
				ps[, ikeep] <- dnorm(graph[, ikeep], mean=m.mat[ikeep,ig], sd=sd.mat[ikeep,ig])
			}
			probs[, ig] <- apply(ps, 1, prod)*priors[ig]
		}
	} else {
		#xbar - sample mean
		#v - sample variance (divided by n not n-1)
		# alpha = xbar(xbar(1-xbar)/v - 1)
		# beta = (1-xbar)(xbar(1-xbar)/v - 1)
		# if in (L, H) xbar = (xbar-L)/(H-L), v = v/(H-L)^2
		# if in (-1, 1) xbar = (xbar +1)/2, v = v/4
		vars <- (v.mat*(nsubj-1)/nsubj)
		xv <- m.mat*(1-m.mat)/vars - 1
		alpha <- m.mat*xv
		beta <- (1-m.mat)*xv
		for (ig in 1:ngroups){
			ps <- matrix(nrow=nsubj, ncol=nkeep)
			for (ikeep in 1:nkeep){
				#ps is the probability given the mean/sd of each group
				ps[, ikeep] <- dbeta(graph[, ikeep], shape1=alpha[ikeep,ig], shape2=beta[ikeep,ig])
			}
			probs[, ig] <- apply(ps, 1, prod)*priors[ig]
		}		
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

ylabels <- rbinom(1000, 5, 0.5)
graph <- agraph <- array(runif(1000 * 100, min=-1, max=1), dim=c(10, 10, 1000))
x <- wssgraph(agraph, ylabels=ylabels, corr=TRUE)
x <- wssgraph(agraph, ylabels=ylabels, corr=FALSE)


graph <- agraph <- array(rnorm(1000 * 100), dim=c(10, 10, 1000))
x <- wssgraph(agraph, ylabels=ylabels)
x <- wssgraph(agraph, ylabels=ylabels, alternative="less")
x <- wssgraph(agraph, ylabels=ylabels, alternative="great")

graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
x <- bssgraph(graph, ylabels=ylabels)
x <- bssgraph(graph, ylabels=ylabels, alternative="less")
x <- bssgraph(graph, ylabels=ylabels, alternative="great")
graph <- data.frame(cbind(graph, ylabels=ylabels))
x <- bssgraph(graph, ylabels=ylabels)
x <- bssgraph(graph, ylabels=ylabels, alternative="less")
x <- bssgraph(graph, ylabels=ylabels, alternative="great")



rm(list=ls())
source("convert.graph.R")
source("get.probs.R")
#### Test data
graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 5, 0.5)
#graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
#graph <- data.frame(cbind(graph, ylabels=ylabels))
upper <- FALSE
alternative = "two.sided"
keep.pct=0.1
weight <- TRUE
corr <- FALSE
workspace=400000

#binary signal subgraph
sig.sgraph <- function(graph, ylabels=NULL, corr=FALSE, upper=FALSE, alternative = "two.sided", keep=0.5, weight=FALSE, workspace=400000, ...){
	
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

	## get levels of adjacencies
	adj <- sort(unique(c(graph)))
	nadj <- length(unique(c(graph)))
	
	if (nadj > 2 & !weight) {
		print("Too many labels - graph is not in {0, 1}, running weighted version")
		weight <- TRUE
	}
	
	if (nadj == 2 & weight) {
		print("Only two labels - using binary classifier")
		weight <- FALSE
	}
	
	## Get number of groups
	ngroups <- length(unique(ylabels))
	if (ngroups < 2) stop("Need at least 2 groups")
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
	pvals <- apply(graph, 2, function(x) getpvals(x, ylabels=ylabels, workspace=workspace, weight=weight))
	
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

	res <- get.probs(graph=graph, weight=weight, corr=corr, ngroups=ngroups, ylabels=ylabels, levs=levs, nkeep=nkeep, nsubj= nsubj)
	probs <- res$probs
	priors <- res$priors
	
	## divide by total to get actual probability and not just numerator P(G = g | Y = y)
	## apply(probs, 1, sum) = P(G=g|Y=1)P(Y=1)+ ... + P(G=g|Y=y_n)P(Y=y_n)
	probs <- probs / apply(probs, 1, sum)
	
	preds <- apply(probs, 1, function(x) which(x == max(x)))
	if (class(preds) != "integer") stop("Problem with Prediction")
	preds <- levs[preds]
	print(table(preds, ylabels))
	print(mean(preds == ylabels))
	
	return(list(ssgraph=sig.graph, priors=priors, prob=probs, pred=preds, ylabels=ylabels))

}


ylabels <- rbinom(1000, 5, 0.5)
ylabels2 <- rbinom(1000, 1, 0.5)

graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
x <- sig.sgraph(agraph, ylabels=ylabels)
x <- sig.sgraph(agraph, ylabels=ylabels2)


graph <- agraph <- array(runif(1000 * 100, min=-1, max=1), dim=c(10, 10, 1000))
x <- sig.sgraph(agraph, ylabels=ylabels, corr=TRUE, weight=TRUE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=TRUE, weight=TRUE)

x <- sig.sgraph(agraph, ylabels=ylabels, corr=TRUE, weight=FALSE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=TRUE, weight=FALSE)

x <- sig.sgraph(agraph, ylabels=ylabels, corr=FALSE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=FALSE)

x <- sig.sgraph(agraph, ylabels=ylabels, corr=FALSE, weight=TRUE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=FALSE, weight=TRUE)

x <- sig.sgraph(graph, ylabels=ylabels, alternative="less")
x <- sig.sgraph(graph, ylabels=ylabels, alternative="great")

graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
graph <- data.frame(cbind(graph, ylabels=ylabels))

x <- sig.sgraph(graph)
x <- sig.sgraph(graph, ylabels=ylabels, alternative="less")
x <- sig.sgraph(graph, ylabels=ylabels, alternative="great")

graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 5, 0.5)
x <- sig.sgraph(agraph, ylabels=ylabels, workspace=10000000)
x <- sig.sgraph(agraph, ylabels=ylabels, alternative="less")
x <- sig.sgraph(agraph, ylabels=ylabels, alternative="great")


ylabels <- c(rep(0, 100), rep(1, 100))
mat0 <- matrix(runif(400, min=0, max=1), nrow=20)
mat1 <- matrix(runif(400, min=0, max=1), nrow=20)

p0 <- c(mat0)
p1 <- c(mat1)

graph0 <- sapply(p0, function(x) rbinom(100, 1, x))
graph1 <- sapply(p1, function(x) rbinom(100, 1, x))

graph <- data.frame(rbind(graph0, graph1))
graph$ylabels <- ylabels
x <- sig.sgraph(graph)




## Rpart 
#fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis)
# unique(rownames(fit$splits))

#randomForest
#iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE, proximity=TRUE)
# iris.rf$importance


#Possible function for formula
# sig.sgraph.formula <- function (formula, data, subset, na.action, ...) 
# {
    # if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
        # "term.labels")) != 1L)) 
        # stop("'formula' missing or incorrect")
    # m <- match.call(expand.dots = FALSE)
    # if (is.matrix(eval(m$data, parent.frame()))) 
        # m$data <- as.data.frame(data)
    # m[[1L]] <- as.name("model.frame")
    # m$... <- NULL
    # mf <- eval(m, parent.frame())
    # DNAME <- paste(names(mf), collapse = " by ")
    # names(mf) <- NULL
    # response <- attr(attr(mf, "terms"), "response")
    # g <- factor(mf[[-response]])
    # if (nlevels(g) != 2L) 
        # stop("grouping factor must have exactly 2 levels")
    # DATA <- split(mf[[response]], g)
    # names(DATA) <- c("x", "y")
    # y <- do.call("t.test", c(DATA, list(...)))
    # y$data.name <- DNAME
    # if (length(y$estimate) == 2L) 
        # names(y$estimate) <- paste("mean in group", levels(g))
    # y
# }


# function (formula, data = NULL, ..., subset, na.action = na.fail) 
# {
    # if (!inherits(formula, "formula")) 
        # stop("method is only for formula objects")
    # m <- match.call(expand = FALSE)
    # if (any(c("xtest", "ytest") %in% names(m))) 
        # stop("xtest/ytest not supported through the formula interface")
    # names(m)[2] <- "formula"
    # if (is.matrix(eval(m$data, parent.frame()))) 
        # m$data <- as.data.frame(data)
    # m$... <- NULL
    # m$na.action <- na.action
    # m[[1]] <- as.name("model.frame")
    # m <- eval(m, parent.frame())
    # y <- model.response(m)
    # Terms <- attr(m, "terms")
    # attr(Terms, "intercept") <- 0
    # m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)), 
        # data.frame(m))
    # for (i in seq(along = ncol(m))) {
        # if (is.ordered(m[[i]])) 
            # m[[i]] <- as.numeric(m[[i]])
    # }
    # ret <- randomForest(m, y, ...)
    # cl <- match.call()
    # cl[[1]] <- as.name("randomForest")
    # ret$call <- cl
    # ret$terms <- Terms
    # if (!is.null(attr(m, "na.action"))) 
        # ret$na.action <- attr(m, "na.action")
    # class(ret) <- c("randomForest.formula", "randomForest")
    # return(ret)
# }


# function (formula, data = NULL, ..., subset, na.action = na.omit, 
    # scale = TRUE) 
# {
    # call <- match.call()
    # if (!inherits(formula, "formula")) 
        # stop("method is only for formula objects")
    # m <- match.call(expand.dots = FALSE)
    # if (identical(class(eval.parent(m$data)), "matrix")) 
        # m$data <- as.data.frame(eval.parent(m$data))
    # m$... <- NULL
    # m$scale <- NULL
    # m[[1]] <- as.name("model.frame")
    # m$na.action <- na.action
    # m <- eval(m, parent.frame())
    # Terms <- attr(m, "terms")
    # attr(Terms, "intercept") <- 0
    # x <- model.matrix(Terms, m)
    # y <- model.extract(m, "response")
    # attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")
    # if (length(scale) == 1) 
        # scale <- rep(scale, ncol(x))
    # if (any(scale)) {
        # remove <- unique(c(which(labels(Terms) %in% names(attr(x, 
            # "contrasts"))), which(!scale)))
        # scale <- !attr(x, "assign") %in% remove
    # }
    # ret <- svm.default(x, y, scale = scale, ..., na.action = na.action)
    # ret$call <- call
    # ret$call[[1]] <- as.name("svm")
    # ret$terms <- Terms
    # if (!is.null(attr(m, "na.action"))) 
        # ret$na.action <- attr(m, "na.action")
    # class(ret) <- c("svm.formula", class(ret))
    # return(ret)
# }

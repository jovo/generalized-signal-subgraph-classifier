rm(list=ls())
source("convert.graph.R")
source("get.probs.R")
source("sig.sgraph.formula.R")
source("sig.sgraph.R")
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


ylabels <- rbinom(1000, 5, 0.5)
ylabels2 <- rbinom(1000, 1, 0.5)
X <- matrix(rnorm(10*1000), nrow=1000)
graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
x <- sig.sgraph(agraph, y=ylabels)
x <- sig.sgraph(agraph, y=ylabels, x=X)
df <- data.frame(cbind(y=ylabels, X))
Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph)
x <- sig.sgraph(agraph, y=ylabels2)


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
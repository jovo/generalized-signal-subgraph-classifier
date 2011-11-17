	
get.probs.formula <- function(graph, x=NULL, y, family=NULL, weight=FALSE, ngroups, levs, nkeep, nsubj, corr=FALSE){
	# if is.null(family)
	
 	eta <- 1/(10*nsubj)
	library(MASS)
	if (!weight & is.null(x)) family <- binomial()
	if (weight & is.null(x)) family <- gaussian()
	# if (weight & is.null(x) & corr=TRUE) family <- beta()
	
	if (weight){
		strfam <- as.character(family$family)
		if (strfam == "gaussian") dens <- "norm"
		if (strfam %in% c("pois", "quasipoisson")) dens <- "pois"
		if (strfam == "binomial") dens <- "binom"
		if (strfam == "Gamma") dens <- "gamma"
		if (strfam == "inverse.gaussian") stop("IG not implemented yet")	
		if (strfam == "multinom") dens <- "multinom"
		if (strfam == "beta") dens <- "beta"
		
		all.probs <- matrix(nrow=nsubj, ncol=ngroups)
		colnames(all.probs) <- levs
		for (ig in 1:ngroups){
			
			## grab the level of ylabels
			ilev <- levs[ig]
			ingroup <- y == ilev
			ingroup[is.na(ingroup)] <- FALSE
			col.probs <- matrix(nrow=nsubj, ncol=ncol(graph))
			for (icol in 1:ncol(graph)){
				run.data <- data.frame(cbind(G=graph[, icol], x))
				### stopped here
				mod <- glm(G ~ ., family=family, data=run.data, subset=ingroup)
				### now using the 
				### Get Ghat  and get parameter estimates from it
				## For Gamma variance/mean = 1/Beta
				if (dens != "multinom") Ghat <- predict(mod, newdata=run.data[ingroup,])
				if (dens == "pois") {
					lambda <- mean(Ghat)
					probs <- dpois(graph, lambda=lambda)
				}
				if (dens == "binom") {
					phat <- mean(Ghat)
					probs <- dbinom(graph, size=1,  prob=phat)
				}
				if (dens == "norm") {
					m <- mean(Ghat)
					sds <- sd(Ghat)
					probs <- dnorm(graph, m, sd=sds)
				}
				if (dens == "gamma") {
					ests <- fitdistr(Ghat, "gamma")
					ests <- ests$estimate
					probs <- dgamma(graph, shape = ests["shape"], rate = ests["rate"])
				}
				if (dens %in% c("beta", "multinom")) {
					stop("Not Implemented yet for this family")
				}
				col.probs[, icol] <- probs

				###need to implement for beta
			}		
			col.probs <- apply(col.probs, 1, prod)
			col.probs <- ifelse(col.probs==1, 1-eta, ifelse(col.probs < .Machine$double.eps, eta, col.probs))
			all.probs[, ig] <- col.probs
		}
	}
return(list(probs=all.probs))
}
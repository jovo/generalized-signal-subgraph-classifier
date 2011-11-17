	
get.probs <- function(graph, weight, corr, ngroups, ylabels, levs, nkeep, nsubj){
	if (weight){
	print("Code chunk 1, weight Probs")	
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
				priors[ig] <- mean(ingroup, na.rm=TRUE)
				## Subset the graphs to only that group
				sset <- graph[ingroup,]
				## get the mean of the adjacency
				## p.mat <- p_{u,v | Y = y}
				m.mat[, ig] <- apply(sset, 2, mean)
				sd.mat[, ig] <- apply(sset, 2, sd)
				v.mat[, ig] <- apply(sset, 2, var)
			}
		print("Code chunk 1 Probs")
		#getprob <- function(graph, priors) {
			# probs 
			probs <- matrix(nrow=nsubj, ncol=ngroups)
			colnames(probs) <- levs
			
			if (corr==FALSE) {
				print("Code chunk 2 Probs")

				for (ig in 1:ngroups){
					ps <- matrix(nrow=nsubj, ncol=nkeep)
					for (ikeep in 1:nkeep){
						#ps is the probability given the mean/sd of each group
						ps[, ikeep] <- dnorm(graph[, ikeep], mean=m.mat[ikeep,ig], sd=sd.mat[ikeep,ig])
					}
					probs[, ig] <- apply(ps, 1, prod)*priors[ig]
				}
			} else {
				print("Code chunk 2 Probs")

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
		} else {

		print("Code chunk 1, Binary Probs")	
	
		##########Get Binary Probabilities############
	
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
			# table
			priors[ig] <- mean(ingroup, na.rm=TRUE)
			## Subset the graphs to only that group
			sset <- graph[ingroup,]
			## get the mean of the adjacency
			## p.mat <- p_{u,v | Y = y}
			p.mat[, ig] <- colMeans(sset)
		}

		print("Code chunk 2, Binary Probs")	
	
		### Can't be probaability = 1 or 0 since dbinom will give 0
		### if htis occurs - may be a good thing - since these may be diagnostic 
		p.mat[p.mat == 1] <- 0.99
		p.mat[p.mat == 0] <- 0.01
	
		# probs 
		probs <- matrix(nrow=nsubj, ncol=ngroups)
		print(levs)
		print(head(probs))
		#colnames(probs) <- levs
		print("Code chunk 3, Binary Probs")	
		
		for (ig in 1:ngroups){
				ps <- log(dbinom(graph, 1, prob=p.mat[,ig]))
				probs[, ig] <- exp(apply(ps, 1, sum))*priors[ig]
		}
		
	}

return(list(probs=probs, priors=priors))
}
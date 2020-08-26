##########################################################################
## Toolkit of usefull function implementations for setting up workspace ##
## when working with WinBUGS from R and for processing MCMC output		##
##########################################################################

###
## A) Setting up workspace

## i) initialize folders which will store MCMC output

initCODAPATH <- function(CODAPATH, BUGSFILENAME) {
	## Delete model .txt file in CODAPATH
	FILE <- paste0(CODAPATH, "/", BUGSFILENAME, ".txt")
	if ( file.exists(FILE) ) {
		file.remove(FILE)
	}
	## Copy the .bugs file to the CODAPATH and rename to .txt extension
	file.copy(from=paste0("./", BUGSFILENAME, ".bug"), to=CODAPATH)
	file.rename(from=paste0(CODAPATH, "/", BUGSFILENAME, ".bug"), to=FILE)
}


###
## B) Loading, processing, and working with MCMC output

## i) Load MCMC output from file

loadMCMC <- function(x, path)
{
	## Loop through all chains and create a list with all the chains
	chains <- list()
	for (i in 1:length(x)) {
		raw <- read.coda(x[i], paste0(path, "/codaIndex.txt"))
		chains[[i]] <- raw
	}
	## Collapse all the chains together
	mcmc <- as.mcmc.list(chains)
	return(as.data.frame(as.matrix(mcmc)))
}

## ii) Functions for computing quantiles of distributions (to use as input to munge function)

## 2.5% quantile
quant025 <- function(x) {
	return(quantile(x, probs=0.025))
}

## 97.5% quantile
quant975 <- function(x) {
	return(quantile(x, probs=0.975))
}


## iii) Functions for plotting MCMC outcome

plotMCMC <- function(prior, posterior, label, outfile, width, height) {
	pdf(file=outfile, width=width, height=height)
	plot(density(prior), type="l", col="blue", lwd=1, lty=2, ylim=c(0, max(c(density(prior)[["y"]], density(posterior)[["y"]]))),
		 xlab="", ylab="", main=paste0("Prior and posterior densities: ", gsub("_", " ", label, fixed=TRUE)), axes=FALSE)
	lines(density(posterior), col="black", lwd=2, lty=1)
	abline(v=0, col="red", lwd=1, lty=2)
	scale <- (10^round(log(as.numeric(q3(c(prior, posterior)) - q1(c(prior, posterior))), 10), 0)) / 2
	lower <- round(as.numeric(quantile(c(prior, posterior), probs=0.001)) / scale, 0) * scale
	upper <- round(as.numeric(quantile(c(prior, posterior), probs=0.999)) / scale, 0) * scale
	axis(1, seq(lower, upper, scale))
	legend("topright", c("Prior","Posterior","Null"), col=c("blue","black","red"), lwd=c(1, 2, 1), lty=c(2, 1, 2))
	graphics.off()
}

## Added an alternate one, due to issue simulating half-cauchy prior
## This one gets the density as input directly

plotMCMC_alt <- function(prior, posterior, label, outfile, width, height) {
	pdf(file=outfile, width=width, height=height)
	plot(prior, type="l", col="blue", lwd=1, lty=2, ylim=c(0, max(c(prior[["y"]], density(posterior)[["y"]]))),
		 xlab="", ylab="", main=paste0("Prior and posterior densities: ", gsub("_", " ", label, fixed=TRUE)), axes=FALSE)
	lines(density(posterior), col="black", lwd=2, lty=1)
	abline(v=0, col="red", lwd=1, lty=2)
	scale <- (10^round(log(as.numeric(q3(c(posterior)) - q1(c(posterior))), 10), 0)) / 2
	lower <- round(as.numeric(quantile(posterior, probs=0.001)) / scale, 0) * scale
	upper <- round(as.numeric(quantile(posterior, probs=0.999)) / scale, 0) * scale
	axis(1, seq(lower, upper, scale))
	legend("topright", c("Prior","Posterior","Null"), col=c("blue","black","red"), lwd=c(1, 2, 1), lty=c(2, 1, 2))
	graphics.off()
}


## iv) Functions for computing Highest Density Intervals (HDI)
## 	   (See Kruschke for rationmale what HDI are better summaries than CrI)

## Sources: http://r.789695.n4.nabble.com/how-to-compute-highest-density-interval-td840831.html
## 			https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
##			https://stackoverflow.com/questions/18635805/how-to-transform-density-object-to-function
## 			(methods: Riemann sum, trapezoidal rule, ecdf. Chose ecdf approach to simplify.)
## (need to test accuracy vs. simulated data; since some steep betas where giving me trouble...)

## Steps: 
## 1: We need a function to compute upper bound based on lower bound
## 	  noting that upper bound will be shifted by a fixed probability relative to lower (e.g., 0.95)
## 2: We need a function to compute the difference between the probability density at the lower
##    and upper bounds. This is used by uniroot to find the point where both densities are the same
##	  which is the first requirement for the HDI
## 3: Finally, we need a function which uses the above two to solve for the lower and upper bounds
##	  of the HDI and return that as the final result
##    Note: for efficiency, we only estimate the ecdf and approximate density functions once

## Deprecated -- use HDIofMCMC, adapted from Kruschke below
getUpperBound <- function(x, Fx, lower, level) {
	## Get cumulative probability mass at lower bound (note: Fx is the ecdf)
	plower <- Fx(lower)
	## Get cumulative probability mass at upper bound
	pupper <- plower + level
	## Compute 'pupper' quantile of the distribution (quantile is the inverse function of ecdf)
	return( as.numeric(quantile(x, probs=pupper)) )
}

## Deprecated -- use HDIofMCMC, adapted from Kruschke below
getDiffDensity <- function(x, Fx, fx, lower, level) {
	## Compute upper bound of HDI (relative to lower bound)
	upper <- getUpperBound(x=x, Fx=Fx, lower=lower, level=level)
	## Return difference in density between lower and upper
	return( abs( fx(lower) - fx(upper) ) )
}

## Deprecated -- use HDIofMCMC, adapted from Kruschke below
getHDI <- function(x, level=0.95, epsilon=1e-9) {
	## Compute empirical cumulative distribution function (ecdf)
	Fx <- ecdf(x)
	## Compute empirical probability density function (epdf)
	fx <- approxfun(density(x), rule=2)
	## Solve for mode of distribution
	fx.mode <- optimize(fx, interval=c(quantile(x, probs=0), quantile(x, probs=1)), maximum=TRUE, tol=1e-9)[["maximum"]]
	## Use uniroot to solve for the lower bound of the HDI
	lowerHDI <- optimize(function(y) { getDiffDensity(x=x, Fx=Fx, fx=fx, lower=y, level=level) },
						 interval=c(quantile(x, probs=0), quantile(x, probs=(1 - level) / 2)), tol=epsilon)[["minimum"]]
	## Retrieve the upper bound of the HDI based on the lower bound we solved for above
	upperHDI <- getUpperBound(x=x, Fx=Fx, lower=lowerHDI, level=level)
	## Return HDI
	return( list(level=level, mode=fx.mode, lowerHDI=lowerHDI, upperHDI=upperHDI) )
}

## Adapted from Kruschke's code from DBDA2
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   a list containing the level, mode, and the limits of the HDI
  
  ## Get HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  
  ## Get mode
  dres = density( sampleVec )
  modeParam = dres$x[which.max(dres$y)]
  
  ## Return results
  return( list(level=credMass, mode=modeParam, lowerHDI=HDImin, upperHDI=HDImax) )
}


## v) Function for plotting nice posterior densities with shaded HDI

plotPosteriorHDI <- function(posterior, label, table, outfile, width, height, nullpt=0) {
	## First, generate kernel density
	fx <- density(posterior)
	## Second, get end point for shaded area (HDI)
	x0 <- table[ table[["REF"]] == label, "lowerHDI"]
	x1 <- table[ table[["REF"]] == label, "upperHDI"]
	## Third, get associated index in density object
	idx0 <- min(which(fx[["x"]] >= x0))
	idx1 <- max(which(fx[["x"]] < x1))
	## Fourth, plot posterior density
	pdf(file=outfile, width=width, height=height)
	plot(fx, type="l", lwd=2,
		 xlab="", ylab="", main=paste0("Posterior density with 95% HDI: ", gsub("_", " ", label, fixed=TRUE)), axes=FALSE)
	## Finally, add shading for HDI
	polygon(x=c(fx[["x"]][c(idx0,idx0:idx1,idx1)]), y=c(0, fx[["y"]][idx0:idx1], 0), col="lightgray")
	abline(v=nullpt, col="red", lwd=1, lty=2)
	scale <- (10^round(log(as.numeric(q3(posterior) - q1(posterior)), 10), 0)) / 2
	lower <- round(as.numeric(quantile(posterior, probs=0.001)) / scale, 0) * scale
	upper <- round(as.numeric(quantile(posterior, probs=0.999)) / scale, 0) * scale
	axis(1, seq(lower, upper, scale))
	graphics.off()
}


## vi) Function for plotting data against posterior predictive quantiles (median, q1, q3, 2.5%, 97.5%)
##     Observations (continuous outcome only) are ordered in ascending order of value with a covariate to second ordering variable

plotPosteriorPredictiveCheck <- function(input, mcmc, y.name, yhat.prefix, byvar, outfile, width, height) {
	## First, order input data by outcome and one chosen covariate
	## And add an id column
	x <- data.frame(id=1:input[["N"]], y=input[[y.name]], x=input[[byvar]])
	x <- x[order(x[["y"]], x[["x"]]), ]
	x[["myindex"]] <- 1:nrow(x)
	## Second, keep chains info needed
	mcmc <- mcmc[ls(mcmc, pattern=yhat.prefix)]
	## Third, prepare the quantiles for each prediction
	stat_table <- as.data.frame(munge(as.data.table(mcmc), fun.list=list(median=median, q1=q1, q3=q3, quant025=quant025, quant975=quant975)))
	stat_table <- stat_table[ grep(yhat.prefix, stat_table[["REF"]]), c("REF","median","q1","q3","quant025","quant975")]
	## Fourth, prepare the plot
	## i) Main plot (observations ordered by outcome, etc.)
	pdf(file=outfile, width=width, height=height)
	plot(x=x[["myindex"]], y=x[["y"]], ylim=c(round(min(stat_table[["quant025"]])), round(max(stat_table[["quant975"]]))),
		 main="Posterior predictive check", xlab="Patient", ylab="Observed outcome")
	## Loop through the 
	for (i in x[["id"]]) {
		## ii) Add in the line segment for the posterior predictive interquartile range
		segments(x[ x[["id"]] == i, "myindex"], stat_table[ stat_table[["REF"]] == paste0(yhat.prefix, i), "q1"],
				 x[ x[["id"]] == i, "myindex"], stat_table[ stat_table[["REF"]] == paste0(yhat.prefix, i), "q3"], col="red")
		## iii) Add in the crosses for the median and 2.5% and 97.5% quantiles
		points(x=rep(x[ x[["id"]] == i, "myindex"], 3), y=stat_table[ stat_table[["REF"]] == paste0(yhat.prefix, i), c("median","quant025","quant975")], col="red", pch=3)
	}
	graphics.off()
}


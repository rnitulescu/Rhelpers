#########################################################################
## Toolkit of usefull function implementations for helping to pin down ##
## prior distributions when working with Bayesian models			   ##
#########################################################################

###
## A) Functions for computing 95% distribution ranges given parameter values
##    (i.e., takes parameters as input and returns a range)

## i) Normal distribution 
getNormalRange <- function(mean, sd) {
	range <- (sd*c(qnorm(0.025), qnorm(0.975))) + mean
	return( list(distribution="normal", lower=range[1], upper=range[2]) )
}

## ii) Uniform

## etc.



###
## B) Functions for pinning down parameter values for various distributions
##	  based on 95% distribution ranges

## i) Normal distribution

## Two equations and two unknowns:
##		z1 * sigma + mu = lower (here z1 is the 2.5 percentile of the standard normal distribution)
##		z2 * sigma + mu = upper (and z2 is the 97.5 percentile of the standard normal distribution)
## We want to solve for mu and sigma, the parameters of the normal distribution that has lower and upper
## as its 2.5 and 97.5 percentiles.
## Thus, in matrix form we get the following:
## ( z1  1 ) ( sigma ) = ( lower )
## ( z2  1 ) ( mu    )   ( upper )
##   2 by 2    2 by 1      2 by 1
##     A         x     =     b
## Solve Ax = b for x (x = solve(a) %*% b)
getNormalParameters <- function(lower, upper, level=0.95) {
	b <- c(lower, upper)
	A <- matrix(c(qnorm((1 - level) / 2), qnorm(1 - ((1 - level) / 2)), 1, 1), 2 ,2)
	x <- round( solve(A) %*% b, 4)
	return( list(mean=x[2], sd=x[1]) )
}

## ii) Log-normal (leverage above implementation)
getLogNormalParameters <- function(lower, upper, level=0.95) {
	b <- log(c(lower, upper))
	params <- getNormalParameters(b[1], b[2], level=level)
	return( list(meanlog=params[["mean"]], sdlog=params[["sd"]]) )
}

## iii) Gamma (tentative: )




###
## C) Plot of priors with 95% CrI
## TODO:


#myroot <- uniroot(function(x) { getPower(r=0.20, RD=0.15, theta=0.5, N=x) - 0.80 }, lower=1e-6, upper=1e+6, tol=1e-9)


#####################################################
## File containing some R function implementations ##
#####################################################

###
## A) Logs

## i) Open and redirect to file
startLog <- function(path, name) {
	FILE <- paste0(path, "/", name, "-", Sys.Date(), ".log")
	if ( file.exists(FILE) ) file.remove(FILE)
	con <- file(FILE)
	sink(con, append=TRUE)
	sink(con, append=TRUE, type="message")
}

## ii) Close and redirect to console
stopLog <- function() {
	sink()
	sink(type="message")
}

###
## B) Data investigation

## i) Print duplicate records (according to byvars)
checkDups <- function(x, byvars) {
	## Collapse data (i.e., count obs)
	x[["nrec"]] <- row.names(x)
	col <- aggregate(as.formula(paste0("nrec", "~", paste(byvars, collapse="+"))), FUN=length, data=x)
	## Retrieve cases where there are 2 or more records per byvars cluster
	col.sub <- col[col[["nrec"]] > 1, ]
	## Merge to main table and return
	x[["nrec"]] <- NULL
	fnl <- merge(x=col.sub[byvars], y=x, by=byvars, all.x=TRUE) ## Left join
	return(fnl)
}

###
## C) Dealing with digits and rounding

## i) Return how many digits we should round to
smart_num_digits <- function(x) {
    x <- abs(x)
    digits <- ifelse(x == 0, 0,
                ifelse(x < 1, ceiling(abs(log(x, 10))),
                  ifelse(x < 10, 2, 
                    ifelse(x < 100, 1, 0))))
    return(digits)
}


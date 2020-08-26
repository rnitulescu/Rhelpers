####################################################################################
## This is the new utils file under development. It will use data.table instead   ##
## of data.frame as the base object. I am refactoring the code into more numerous ##
## and smaller functions and using more efficient algorithms (i.e., lapply, Map,  ##
## etc., instead of 'for' loops)												  ##
####################################################################################

###
## A) We need to require that data.table be loadable (and loaded)
if (!suppressWarnings(require(data.table, quietly=TRUE))) stop("Can't load data.table library.")
if (!suppressWarnings(require(xtable, quietly=TRUE))) stop("Can't load xtable library.")

###
## B) We implement the lower-order functions which will be called by the higher-order ones

## Symbolic constant for list fo numeric variable types
NUMLST <- c("integer","numeric","double")
FACLST <- c("logical","factor")

## Return the minimum value in a vector, excluding NA's and returning NA if all values in vector ar NA's
safe.min <- function(x) {
	tryCatch(min(x, na.rm=TRUE), warning=function(w) {return(NA)})
}

## Return the maximum value in a vector, excluding NA's and returning NA if all values in vector ar NA's
safe.max <- function(x) {
	tryCatch(max(x, na.rm=TRUE), warning=function(w) {return(NA)})
}

## Return the first quartile of the distribution in a vector
q1 <- function(x) {
	if ( !(class(x) %in% NUMLST) ) {
		stop("Invalid variable type")
	} else {
		return(quantile(x, probs=0.25, na.rm=TRUE))
	}
}

## Return the second quartile (i.e., median) of the distribution in a vector
q2 <- function(x) {
	if ( !(class(x) %in% NUMLST) ) {
		stop("Invalid variable type")
	} else {
		return(quantile(x, probs=0.50, na.rm=TRUE))
	}
}

## Return the third quartile of the distribution in a vector
q3 <- function(x) {
	if ( !(class(x) %in% NUMLST) ) {
		stop("Invalid variable type")
	} else {
		return(quantile(x, probs=0.75, na.rm=TRUE))
	}
}

## Return the number of values in a vector
total <- function(x) {
	return(length(x))
}

## Return the number of values in a vector which are not NA's
notna <- function(x) {
	return(length(x[!is.na(x)]))
}

## Return the number of values in a vector which are NA's
na <- function(x) {
	return(length(x[is.na(x)]))
}

## Return a one-way table counting the number of observations per category from a vector
## (logical: TRUE/FALSE, factor: levels) Also computes percentage
tab1way <- function(x) {
	## Don't run if variables are not of logical or factor type
	if ( !(class(x) %in% FACLST) ) {
		stop("Invalid variable type")
	} else {
		## Compute frequency and proportion (but drop the FALSE entries as they are redundant for for logical types)
		if ( class(x) == "logical" ) {
			freq <- as.data.table(table(x, exclude="useNA"))[x != "FALSE"]
			perc <- as.data.table(prop.table(table(x, exclude="useNA")))[x != "FALSE"]
		} else {
			freq <- as.data.table(table(x, exclude="useNA"))
			perc <- as.data.table(prop.table(table(x, exclude="useNA")))
		}
		## Rename columns so we can merge freq and perc
		setnames(freq, old=c("x","N"), new=c("CAT","freq"))
		setnames(perc, old=c("x","N"), new=c("CAT","prop"))
		## Merge freq and perc on CAT
		y <- freq[perc, on="CAT"]
		## Convert proportion to percentage and drop prop column (forcing use of standard evaluation)
		set(y, j="perc", value = 100 * y[["prop"]])
		set(y, j="prop", value = NULL)
		## Return table
		return(y)
	}
}

## Return a two-way table counting the number of observations per category combination from two vectors
## (logical: TRUE/FALSE, factor: levels) Also computes percentage
tab2way <- function(x, byvar) {
	## Don't run if variables are not of logical or factor type
	if ( !(class(x) %in% FACLST) | !(class(byvar) %in% FACLST) ) {
		stop("Invalid variable type")
	## Don't run if there are missing values for the byvar, since this leads to warning in reshape
	} else if ( length(byvar[is.na(byvar)]) > 0 ) {
		stop("<NA> values in byvar")
	} else {
		## Compute frequency and proportion
		freq <- as.data.table(table(x, byvar, exclude="useNA"))
		perc <- as.data.table(prop.table(table(x, byvar, exclude="useNA"), 2))
		## Rename columns so we can merge freq and perc
		setnames(freq, old=c("x","N"), new=c("CAT","freq"))
		setnames(perc, old=c("x","N"), new=c("CAT","prop"))
		## Merge freq and perc on CAT and byvar
		y <- freq[perc, on=.(CAT, byvar)]
		## Convert proportion to percentage and drop prop column (forcing use of standard evaluation)
		set(y, j="perc", value = 100 * y[["prop"]])
		set(y, j="prop", value = NULL)
		## Restructure to have columns for each byvar level
		yt <- reshape(y, v.names=c("freq","perc"), timevar="byvar", idvar="CAT", direction="wide", sep="_")
		return(yt)
	}
}

## Return the arithmetic mean of a numeric vector
armean <- function(x) {
	if ( !(class(x) %in% NUMLST) ) {
		stop("Invalid variable type")
	} else {
		return(mean(x, na.rm=TRUE))
	}
}

## Return the standard deviation of a numeric vector
stdev <- function(x) {
	if ( !(class(x) %in% NUMLST) ) {
		stop("Invalid variable type")
	} else {
		return(sd(x, na.rm=TRUE))
	}
}

###
## C) We implement the higher-order function which will be called by the highest-order ones

## Core function for summarizing columns of data
## Returns a data table object with the summarized data in a standardized format
fn.summary <- function(x, byvar=NULL, cols, FUN, FUN.name)
# x: data table object containing columns to be summarized
# byvar: character variable representing name of variable to group statistics by (default: NULL)
# cols: character vector representing names of variables to compute statistics for
# FUN: function to use for computing statistics
# FUN.name: name of function to use for computing statistics (this function will only be called by another in a well controlled manner, so this will not be an issue)
{
	## Don't run if there are multiple byvars. Implementaton only works with 1 byvar at the moment.
	if ( length(byvar) > 1 ) {
		stop("Invalid arguments")
	} else if ( !(FUN.name %in% c("tab1way","tab2way")) ) {
	## Case 1: Numeric variables (i.e., not for 1-way or 2-way tables)
		if ( is.null(byvar) ) {
			## Summarize data
			y <- x[, lapply(.SD, FUN), .SDcols=cols]
			## Transpose results
			yt <- transpose(y)
			## Rename results column
			setnames(yt, FUN.name)
			## Merge in the names of the variables summarized
			yt <- cbind(data.table(REF=colnames(y)), yt)
		} else {
			## Summarize data
			y <- x[, lapply(.SD, FUN), by=byvar, .SDcols=cols]
			## Transpose results
			yt <- transpose(y)[-1]
			## Format header/column names
			header <- y[[byvar]]
			header <- paste(FUN.name, header, sep="_")
			setnames(yt, header)
			## Merge in the names of the variables summarized and the byvar
			yt <- cbind(data.table(REF=colnames(y)[-1]), byvar, yt)
		}
	} else {
	## Case 2: Deal with the case of categorical (i.e., factor) variables
		if ( is.null(byvar) ) {
		## One-way table
			## Summarize data
			y <- lapply(x[, cols, with=FALSE], FUN)
			## Merge in name of variable
			y <- Map(function(y, i) cbind(VAR=i, y), y, names(y))
		} else {
		## Two-way table
			## Summarize data
			y <- Map(function(z) FUN(z, x[[byvar]]), x[, cols, with=FALSE])
			## Merge in name of variable and byvar
			y <- Map(function(y, i) cbind(VAR=i, byvar, y), y, names(y))
		}
		## Stack tables
		yt <- do.call("rbind", y)
		## Create reference column, drop unneeded columns, and reorder columns
		set(yt, j="REF", value = paste(yt[["VAR"]], yt[["CAT"]], sep="_"))
		set(yt, j=c("VAR","CAT"), value = NULL)
		setcolorder(yt, c("REF", setdiff(names(yt), "REF")))
	}
	return(yt)
}

###
## D) We implement the highest-order functions which end the stack

## Returns class of a column inside a data table
getClass <- function(name, x) {
	return(class(x[[name]])[1])
}

## One-way or two-way analysis of data
## Calls function 'fn.summary' for statistics and stacks all results into one data table which the function returns
## Takes data table as input and byvar variable name. Also, it allows for using a different set of functions than the default, if passed appropriately
munge <- function(x, byvar=NULL, fun.list=list(armean=armean, stdev=stdev, q1=q1, q2=q2, q3=q3, safe.min=safe.min, safe.max=safe.max, total=total, na=na, notna=notna)) {
	## Don't run if there are multiple byvars. Implementaton only works with 1 byvar at the moment.
	if ( length(byvar) > 1 ) {
		stop("Invalid arguments")
	} else {
		## First, group variables by class (excluding the byvar)
		nm.lst <- setdiff(names(x), byvar)
		cl.lst <- unlist(lapply(nm.lst, getClass, x))
		cont.vars <- nm.lst[cl.lst %in% NUMLST]
		disc.vars <- nm.lst[cl.lst %in% FACLST]
		## Second, compute stats (1-way or 2-way, as the case may be)
		## Case 1: 1-way stats
		if ( is.null(byvar) ) {
			## A) compute 1-way frequency tables, if the list of such variables is non-empty
			if ( length(disc.vars) > 0 ) {
				disc.stat <- fn.summary(x, cols=disc.vars, FUN=tab1way, FUN.name="tab1way")
			}
			## B) compute statistics for continuous variables and merge them together, if the list of such variables is non-empty
			if ( length(cont.vars) > 0 ) {
				cont.stat <- Map(function(y, z) fn.summary(x, cols=cont.vars, FUN=y, FUN.name=z), fun.list, names(fun.list))
			}
			## C) Add a row for the number of obs in table
			tmp <- data.table(REF="NOBS", freq=nrow(x))
		} else {
		## Case 2: 2-way stats
			## A) compute 2-way frequency tables, if the list of such variables is non-empty
			if ( length(disc.vars) > 0 ) {
				disc.stat <- fn.summary(x, byvar=byvar, cols=disc.vars, FUN=tab2way, FUN.name="tab2way")
			}
			## B) compute statistics for continuous variables and merge them together, if the list of such variables is non-empty
			if ( length(cont.vars) > 0 ) {
				cont.stat <- Map(function(y, z) fn.summary(x, byvar=byvar, cols=cont.vars, FUN=y, FUN.name=z), fun.list, names(fun.list))
			}
			## C) Add a row for the number of obs in table
			## Tabulate number of observations per group (as per byvar)
			tmp <- transpose(x[, .N, by=byvar])
			## Define column names based on contens of first row of this table (which contains the names of the byvar groups)
			header <- paste0("freq_", as.matrix(tmp[1]))
			## Set the names of the columns
			setnames(tmp, header)
			## Keep on the second row (i.e., the row containing the counts)
			tmp <- tmp[2]
			## Add the structural columns (REF and byvar)
			tmp <- cbind(REF="NOBS", byvar=byvar, tmp)
			## Change values back to numeric (since they were forced to character when we transposed above)
			tmp <- tmp[, lapply(.SD, as.numeric), by=c("REF","byvar")]
		}
		## Third, merge together all tables in list of continuous variables, if the list of such variables in non-empty
		if ( length(cont.vars) > 0 ) {
			cont.stat.mrg <- Reduce(merge, cont.stat)
			## Change values back to numeric (since they were forced to character when we Mapped the summary functions above)
			if ( is.null(byvar) ) {
				cont.stat.mrg <- cont.stat.mrg[, lapply(.SD, as.numeric), by="REF"]
			} else {
				cont.stat.mrg <- cont.stat.mrg[, lapply(.SD, as.numeric), by=c("REF","byvar")]
			}
		}
		## Finally, return stacked stats table, depending on what was computed
		if ( length(disc.vars) == 0 & length(cont.vars) == 0 ) {
			return(tmp)
		} else if ( length(disc.vars) == 0 ) {
			return( rbindlist(list(tmp, cont.stat.mrg), fill=TRUE) )
		} else if ( length(cont.vars) == 0 ) {
			return( rbindlist(list(tmp, disc.stat), fill=TRUE) )
		} else {
			return( rbindlist(list(tmp, disc.stat, cont.stat.mrg), fill=TRUE) )
		}	
	}
}


## Finally, we format results, whether one-way or two-way (i.e., unstratified or stratified)
## Function takes the data table containing munged output as well as an optional affix (stratified output only)
## Returns data frame with one column representing the formatted results (default is now just freq (perc%) or median (Q1 ; Q3))
formatMungedData <- function(x, affix="") {
	## Convert data table to data frame
	x <- as.data.frame(x)
	## Set up names of columns (with or without an affix) for both input and output
	## (i.e., when data was munged stratified we need to track statum using an affix)
	if (affix == "") {
		column.name <- "frpc_med_iqr"
		sfreq <- "freq"; sperc <- "perc"
		sq1 <- "q1"; sq2 <- "q2"; sq3 <- "q3"
	} else {
		column.name <- paste("frpc_med_iqr", affix, sep="_")
		sfreq <- paste("freq", affix, sep="_"); sperc <- paste("perc", affix, sep="_")
		sq1 <- paste("q1", affix, sep="_");	sq2 <- paste("q2", affix, sep="_")
		sq3 <- paste("q3", affix, sep="_")
	}
	## Create column (format data for presentation)
	x[[column.name]] <- ifelse(!is.na(x[[sfreq]]),
						   ifelse(!is.na(x[[sperc]]), sprintf("%.0f (%.0f%%)", x[[sfreq]], x[[sperc]]), x[[sfreq]]),
						   ifelse(!is.na(x[[sq2]]), sprintf("%.0f (%.0f ; %.0f)", x[[sq2]], x[[sq1]], x[[sq3]]), NA))
	## Return column of interest
	return(x[column.name])
}


###
## E) Table exporting

## i) For one-way analyses, this exports to latex
munge_to_latex <- function(x, prefix=NULL, mycaption="", outpath=NULL, tab.envir="tabular", regex.mode=1, float.option=TRUE) {
    ## Step 0: Check input
    if (is.null(prefix) & is.null(outpath)) stop("Please define prefix and outfile")
    if (!(regex.mode %in% c(1, 2))) stop("Invalid regex.mode")
    if (!(float.option %in% c(TRUE, FALSE))) stop("Invalid float.option")
    ## And create outpath if it does not already exist
    if (!dir.exists(outpath)) dir.create(outpath)

    ## Save input for later use and reuse
    input <- x

    ## Step 1: Start with categorical variables

    ## If we can (i.e., if "perc" in the table), then do
    if (any(names(input) %in% "perc")) {
        ## i) Extract columns we care about    
        x <- as.data.frame(input)
        x <- x[ !is.na(x[["freq"]]) , c("REF","freq","perc") ]
        ## ii) Clean up reference column
        x[ x[["REF"]] == "NOBS" , "REF" ] <- "Number of observations"
        ## iii) Format results
        x[["n (%)"]] <- ifelse(!is.na(x[["perc"]]), 
                               sprintf("%d (%.0f%%)", x[["freq"]], x[["perc"]]),
                               sprintf("%d", x[["freq"]]))
        ## iv) Start preparing final table (below processing will be to add formating)
        y <- x
        if (regex.mode == 1) {
            y[["Variable"]] <- gsub("^.*_([[:alnum:]]+)$", "\\1", trimws(y[["REF"]]))
        } else {
            y[["Variable"]] <- gsub("^.*_(.+)$", "\\1", trimws(y[["REF"]]))
        }
        y <- y[, c("Variable","n (%)") ]
        ## v) Get rid of underscores in variable names
        y[["Variable"]] <- gsub("_", " ", y[["Variable"]], fixed=TRUE)

        ## vi) Truncate group label from categorical variables
        if (regex.mode == 1) {
            x[["REF2"]] <- gsub("^(.*)_[[:alnum:]]+$", "\\1", trimws(x[["REF"]]))
        } else {
            x[["REF2"]] <- gsub("^(.*)_.+$", "\\1", trimws(x[["REF"]]))
        }
        ## vii) Identify where we will insert additional formating for groups
        to_insert <- as.list(as.numeric(row.names(unique(x[c("REF2")]))[-1])-1)
        ## viii) Define the addtorow list for passing to print method of xtable
        addtorow <- list(pos=to_insert, 
                         command=paste0("\\hline\n",
                                        "\\emph{", gsub("_", " ", unique(x[["REF2"]])[-1], fixed=TRUE), "} & \\\\\n",
                                        "\\hline\n"))

        ## ix) Finally, we can now print the final table with the additional formating
        print(xtable(y, align="|l|l|c|", type="latex",
                     caption=paste(mycaption, ": categorical variables")),
              tabular.environment=tab.envir,
              include.rownames=FALSE,
              add.to.row=addtorow, floating=float.option,
              file=paste(outpath, paste(prefix, "categorical.tex", sep="_"), sep="/"))
    }


    ## Step 2: Now deal with continuous variables

    ## If we can (i.e., if "armean" in the table), then do
    if (any(names(input) %in% "armean")) {
        ## i) Extract columns we care about
        x <- as.data.frame(input)
        x <- x[ is.na(x[["freq"]]) , setdiff(names(x), c("freq","perc")) ]
        ## ii) Format results
        x[["Mean (SD)"]] <- ifelse(!is.na(x[["stdev"]]), 
                                   sprintf(paste0("%.", smart_num_digits(x[["armean"]]), 
                                                  "f (%.", smart_num_digits(x[["stdev"]]), "f)"),
                                           x[["armean"]], x[["stdev"]]),
                                   sprintf(paste0("%.", smart_num_digits(x[["armean"]]), "f"), x[["armean"]]))
        x[["Median (Q1, Q3)"]] <- ifelse(!is.na(x[["q1"]]), 
                                   sprintf(paste0("%.", smart_num_digits(x[["q2"]]), 
                                                  "f (%.", smart_num_digits(x[["q1"]]), "f, ",
                                                  "%.", smart_num_digits(x[["q3"]]), "f)"),
                                           x[["q2"]], x[["q1"]], x[["q3"]]),
                                   sprintf(paste0("%.", smart_num_digits(x[["q2"]]), "f"), x[["q2"]]))
        x[["Min"]] <- sprintf(paste0("%.", smart_num_digits(x[["safe.min"]]), "f"), x[["safe.min"]])
        x[["Max"]] <- sprintf(paste0("%.", smart_num_digits(x[["safe.max"]]), "f"), x[["safe.max"]])
        x[["NA's"]] <- sprintf(paste0("%.0f"), x[["na"]])

        ## iii) Keep what we want to export
        names(x)[names(x) == "REF"] <- "Variable"
        x <- x[c("Variable","Mean (SD)","Median (Q1, Q3)","Min","Max","NA's")]
        
        ## iv) Get rid of underscores in variable names
        x[["Variable"]] <- gsub("_", " ", x[["Variable"]], fixed=TRUE)

        ## v) Finally, we can now print the final table
        print(xtable(x, align="|l|l|ccccc|", type="latex",
                     paste(mycaption, ": continuous variables")),
              tabular.environment=tab.envir,
              include.rownames=FALSE, floating=float.option,
              file=paste(outpath, paste(prefix, "continuous.tex", sep="_"), sep="/"))
    }
}


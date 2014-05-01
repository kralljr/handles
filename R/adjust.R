#' Adjusting data below the MDL
#'
#' \code{adjust} adjusts censored concentrations below the MDL
#'
#' These are functions to create complete PM2.5 constitutent
#' concentrations by substituting, excluding, or imputing censored 
#' data. Works on dataframe where the first column is date and all subsequent 
#' columns are concentrations of chemical constituents.
#'
#' @title adjust
#' @param data data frame of daily constituent concentrations with date as first column
#' @param mdl vector or matrix of MDLs for each variable of data
#' @param adjust one of three methods.  'substitute' substitutes all censored constant proportion of the mdl (sub * mdl).  'exclude' excludes constituents with more than experc censored data.  'likelihood' multiply imputes censored data with a likelihood-based method.
#' @param sub proportion of MDL to substitute censored concentrations.  Default is 0.5.
#' @param experc If a constituent exceeds this value, using the exclude method, the constituent will be dropped from the analysis. 
#' @param N number of draws from posterior for likelihood-based method
#' @param burnin number of samples to discard for likelihood-based method
#' #param ... other arguments
#' @export
#' @examples
#' data(nycdat)
#' data(nycmdl)
#' #fix data totals
#' pm <- nycdat$PM25
#' whPM <- which(colnames(nycdat) == "PM25")
#' nycdat <- nycdat[, -whPM]
#' whPM <- which(colnames(nycmdl) == "PM25")
#' nycmdl <- nycmdl[, -whPM]
#' adjust(nycdat, mdl = nycmdl, adjust = "substitute", sub = 1/sqrt(2))
#' adjust(nycdat, mdl = nycmdl, adjust = "exclude", experc = 0.2)
#' adjust(nycdat, mdl = nycmdl, adjust = "likelihood")
adjust <- function(x, ...) UseMethod("adjust")


#' @rdname adjust
#' @export
adjust.default <- function(data, mdl, 
	adjust = c("substitute", "exclude", "likelihood"),
	sub = 0.5, experc = 0.25, N = 10, burnin = 2, ...) {
		
	if(class(data[, 1]) != "Date") {
		stop("First column must be 'date'")
	}
	dates <- data[, 1]
	dat <- as.matrix(data[, -1])
	cn <- colnames(dat)
	
	
	#get one MDL for each constituent	
	if(!is.null(dim(mdl)) & dim(mdl)[1] > 1) {
		var1 <- apply(mdl, 2, var, na.rm = T)
		if(sum(var1) > 0) {
			warning("MDLs vary over time, using maximum")
		}
		mdl <- apply(mdl, 2, max, na.rm = T)
	}
	mdl <- matrix(rep(as.numeric(mdl), nrow(dat)), 
		byrow = T, nrow(dat))
	
	bdls <- 1 * (dat < mdl)
	sumbdls <- apply(bdls, 2, mean)
	
	out <- list()
	out$call <- match.call()
	out$bdls <- sumbdls
	
	if(adjust == "substitute") {
		dat1 <- dat * (1- bdls) + sub * mdl * bdls
		out$sub <- sub
	}else if(adjust == "exclude") {
		
		whkeep <- which(sumbdls <= experc)
		if(length(whkeep) == 0) {
			stop("Exclude method drops all constituents")
		}
		
		dat1 <- dat[, whkeep]
		out$exclude <- cn[-whkeep]
		cn <- cn[whkeep]
	}else if(adjust == "likelihood") {
		guess <- list()
		guess[[1]] <- dat
		guess[[2]] <- colMeans(dat)
		guess[[3]] <- cov(dat)
		lhood1 <- lhood(dat, mdl, guess, burnin, N)
		out$impmean <- lhood1$gthet
		out$impcov <- lhood1$gsig
		out$impdat <- lhood1$gymiss
		out$burnin <- burnin
		out$N <- N
		dat1 <- apply(lhood1$gymiss, c(1, 2), mean, na.rm = T)
	}
	
	dat1 <- data.frame(dates, dat1)
	colnames(dat1) <- c("Date", cn)
	out$dat <- dat1
	out$adjust <- adjust
	
	class(out) <- "adjust"
	out
}



#' @rdname adjust
#' @export
print.adjust <- function(x) {
	cat("Call:\n")
	print(x$call)
	cat("Head of adjusted data:\n")
	print(head(x$dat))
	cat("Proportion below the MDL:\n")
	print(x$bdls)
	
	if(x$adjust == "substitute") {
		cat("Substitute below MDL with MDL *\n")
		print(x$sub)
	}else if(x$adjust == "exclude") {
		cat("Constituents excluded:\n")
		print(x$exclude)
	}else if(x$adjust == "likelihood") {
		cat("Samples: N =",  x$N, ", Burnin:", x$burnin, "\n")
	}
}
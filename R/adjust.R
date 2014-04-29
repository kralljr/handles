


adjust <- function(x, ...) UseMethod("adjust")


adjust.default <- function(data, mdl, 
	method = c("substitute", "exclude", "likelihood"),
	sub = 0.5, experc = 0.25, N = 10, burnin = 2, ...) {
		
	if(class(data[, 1]) != "Date") {
		stop("First column must be 'date'")
	}
	dates <- data[, 1]
	dat <- data[, -1]
	
	#get one MDL for each constituent	
	if(!is.null(dim(mdl))) {
		var1 <- apply(mdl, 2, var, na.rm = T)
		if(sum(var1) > 0) {
			warning("MDLs vary over time, using maximum")
		}
		mdl <- apply(mdl, 2, max, na.rm = T)
	}
	mdl <- matrix(rep(mdl, nrow(dat)), byrow = T, nrow(dat))
	
	bdls <- 1 * (dat < mdl)
	sumbdls <- apply(bdls, 2, mean)
	
	out <- list()
	out$call <- match.call()
	out$bdls <- sumbdls
	
	if(method == "substitute") {
		dat1 <- dat * (1- bdls) + sub * mdl * bdls
		out$sub <- sub
	}else if(method == "exclude") {
		cn <- colnames(dat)
		whkeep <- cn[which(sumbdls <= experc)]
		if(length(whkeep) == 0) {
			stop("Exclude method drops all constituents")
		}
		
		dat1 <- subset(dat, subset = cn %in% whkeep)
		out$exclude <- subset(cn, subset = !(cn %in% whkeep))
	}else if(method == "likelihood") {
		guess <- list()
		guess[[1]] <- dat
		guess[[2]] <- colMeans(dat)
		guess[[3]] <- cov(dat)
		adjust <- lhood(dat, mdl, guess, burnin, N)
		out$impmean <- adjust$gthet
		out$impcov <- adjust$gsig
		out$impdat <- adjust$gymiss
		out$burnin <- burnin
		out$N <- N
		dat1 <- apply(adjust$gymiss, c(1, 2), mean, na.rm = T)
	}
	
	dat1 <- data.frame(dates, dat1)
	colnames(dat1) <- c("Date", colnames(dat))
	out$dat <- dat1
	out$method <- method
	
	class(out) <- "adjust"
	out
}



print.adjust <- function(x, ...) {
	cat("Call:\n")
	print(x$call)
	cat("Head of adjusted data:\n")
	print(head(x$dat))
	cat("Proportion below the MDL:\n")
	print(x$bdls)
	
	if(x$method == "substitute") {
		cat("Substitute below MDL with MDL *\n")
		print(x$sub)
	}else if(x$method == "exclude") {
		cat("Constituents excluded:\n")
		print(x$exclude)
	}else if(x$method == "likelihood") {
		cat("Samples: N =",  x$N, ", Burnin:", x$burnin, "\n")
	}
}
#' Absolute Principal Component Analysis (APCA)
#'
#' \code{simAPCA} applies APCA and extracts summary information for simulation study
#'
#' @param data data frame of daily constituent concentrations with date as first column
#' @param tots vector of total concentrations (total PM2.5) for each day.  If null, uses \code{rowSums(data)}
#' @param nsources number of sources.  If null, uses number of eigenvalues of the correlation matrix greater than one.
#' @param adjust method to adjust censored concentrations below minimum detection limits (MDLs).  If null, uses complete case data
#' @param mdl either a vector of mdls corresponding to each constituent or a matrix of mdls corresponding to each constituent and each day. 
#' @param cut cutoff for eigenvalues.  Only used if nsources is NULL
#' @param type Source apportionment method (apca or mapca)
#' @param mons vector of monitors corresponding to data
#' @param i monitor number (for printing)
#' @param print whether to print monitors without temporal variability
#' @export
SIMapca <- function(data, tots = NULL,
	nsources = NULL, adjust = NULL,  
	mdl = NULL, cut = 1, type = "apca", mons = NULL,
    i = NA, print = F, k = 10, complete = F,...) {
	
	
	apca1 <- apca(data, tots, nsources, adjust, mdl, cut,
		type, mons, i, print, ...)
		
	#summary
	conc <- apca1$conc
	means <- apply(conc, 2, mean, na.rm = T)	
	sd <- apply(conc, 2, sd, na.rm = T)	
		
	#classify
	profile <- apca1$prof
	data(speciate)
	class <- classify(profile = profile, k = k, 
		sprofiles = speciate)
		
	out <- list(means = means, sd = sd, class = class)	
	if(complete) {
		out[["apca"]] <- apca1
	}
	
	out
		
	
}



#function to adjust sprofiles to include subset
# of cons
# sprofiles: restricted speciate database
# profile: profile of data
dropcons <- function(sprofiles, profile) {
	sprofiles <- sprofiles[, colnames(profile)]
	
	#rescale
	rsums <- rowSums(sprofiles)
	
	#if divide by zero issue
	wh0 <- which(rsums == 0)
	if(length(wh0) > 0) {
		sprofiles <- sprofiles[-wh0, ]
		rsums <- rsums[-wh0]
	}
	
	sprofiles <- sweep(sprofiles, 1, rsums, "/")
	sprofiles
	
	
}





#' \code{classify} classifies profiles according to speciate for simulation study
#'
#' @param profile Source apportionment profile matrix with number of columns equal to number of chemical constituents and number of rows equal to number of sources
#' @param k k corresponding to k nearest neighbors
#' @param sprofiles database of profiles from speciate
#' @export
classify <- function(profile, k = 10, 
	sprofiles = speciate){

	nsource <- nrow(profile)
	
	#if we have dropped some constituents
	if(ncol(sprofiles) != ncol(profile)) {
		sprofiles <- dropcons(sprofiles, profile)
	}
	
	
	#get labels
	labels <- factor(rownames(sprofiles))
	
	#classify according to k nearest neighbors
	class <- knn(sprofiles, profile, labels,
		k = k, prob = T)
	probs <- attributes(class)$prob
	class <- as.character(class)
	
	
	#remove duplicate classification
	class <- rmdups(profile = profile, class = class, 
		probs = probs, k = k, sprofiles = sprofiles)
	
	class	

}





#' \code{rmdups} reassigns classification for duplicate classifications
#'
#' @param profile Source apportionment profile matrix with number of columns equal to number of chemical constituents and number of rows equal to number of sources
#' @param class character vector of classifications
#' @param probs vector of probabilities corresponding to classifications
rmdups <- function(profile, class, probs, k,
	sprofiles = speciate) {
	
	#number of sources
	nsource <- nrow(profile)

	#name class to preserve original order
	seqs <- seq(1, nsource)
	names(class) <- seqs
	
	#reorder class/probs/profile by probs
	class <- class[order(probs, decreasing = T)]
	probs <- probs[order(probs, decreasing = T)]
	profile <- profile[order(probs, decreasing = T), ]


	#while duplicated sources
	while(length(unique(class)) != nsource) {

		#find first duplicated
		mns <- min(which(duplicated(class)))
		
		#remove classes already classified
		keeps <- class[1 : (mns - 1)]
		rn <- factor(rownames(sprofiles))
		sprofiles <- sprofiles[-which(rn %in%  keeps), ]
		rn <- factor(rownames(sprofiles))

		#which need to be reclass
		reclass <- c(mns : nsource)
		#reclassify remaining
		test <- knn(sprofiles, profile[reclass, ], 
			cl = rn, k = k, prob = T)
		probs[reclass] <- attributes(test)$prob
		class[reclass] <- as.character(test)


		#reorder reclass by probability
		classt <- class[reclass][order(probs[reclass], 
			decreasing = T)]
		class[reclass] <- classt
		
		
		names(class) <- c(names(keeps), names(classt))

		}

	#reorder by original names
	class <- class[order(names(class))]
	
	class
	
}


		




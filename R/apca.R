#' Absolute Principal Component Analysis (APCA)
#'
#' \code{apca} performs APCA on daily PM2.5 constituent concentrations
#'
#' This is a function to estimate PM2.5 source profiles and daily 
#' PM2.5 source concentrations from PM2.5 constituent concentrations
#' observed at ambient monitors (e.g. EPA Chemical Speciation Network).
#' Works on dataframe where the first column is date and all subsequent 
#' columns are concentrations of chemical constituents.
#' See Thurston and Spengler (1985, Atmospheric Environment)
#'
#' @param data data frame of daily constituent concentrations with date as first column
#' @param tots vector of total concentrations (total PM2.5) for each day.  If null, uses \code{rowSums(data)}
#' @param nsources number of sources.  If null, uses number of eigenvalues of the correlation matrix greater than one.
#' @param adjust method to adjust censored concentrations below minimum detection limits (MDLs).  If null, uses complete case data
#' @param mdl either a vector of mdls corresponding to each constituent or a matrix of mdls corresponding to each constituent and each day. 
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
#' apca(nycdat)
#' apca(nycdat, tots = pm)
#' apca(nycdat, mdl = nycmdl, adjust = "substitute")
apca <- function(x, ...) UseMethod("apca")

#' @rdname apca
#' @export
apca.default <- function(data, tots = NULL,
	nsources = NULL, adjust = NULL,  
	mdl = NULL, cut = 1, type = "apca", mons = NULL,
    i = NA, print = F, ...){
		
	if(!is.null(adjust)) {
		adj1 <- adjust(data = data, mdl = mdl, 
			adjust = adjust, ...)
		data <- adj1$dat
	}	
		
	if(class(data[, 1]) != "Date") {
		stop("First column must be 'date'")
	}
		
	dates <- data[, 1]
	data <- data[, -1]	
		
	#standardize data (mean zero, sd 1)
	stddata <- stdize1(data, i, print)[[1]]
	
	#Use PC method to get number of sources
	if(is.null(nsources)) {
		pr1 <- prcomp(stddata)
		nsources <- length(which(pr1$sdev > cut))
	}
		
	#get scores and rotation matrix bstar	
	temp <- getbstar(stddata = stddata, nsources = nsources, 
		dates = dates)
	scores <- temp$scores
	bstar <- temp$bstar
	vmax <- temp$vmax

	#use regressions to obtain source profiles and source concentrations
    apca <- getscoresprofs(data = data, bstar = bstar, 
	    	scores = scores, tots = tots, dates = dates, mons = mons, type = type) 
    
	apca$vmax <- vmax
	apca$dates <- dates
	apca$nsources <- nsources
	apca$call <- match.call()
	

	if(!is.null(adjust)) {
		apca$adjust <- adj1
	}	

	class(apca) <- "apca"
	apca
}





	
	
#' \code{stdize1} Standardizes matrix so each column is mean 0 with variance 1
#'
#' @param data data frame of daily constituent concentrations with date as first column
#' @param i iteration number
#' @param print whether to print if variability of column is 0
#' @export
stdize1 <- function(data, i = "NA", print = F) {
	
	#get means and sd
	cm <- colMeans(data, na.rm = T)
	sds <- apply(data, 2, sd, na.rm = T)
	
	#standardize
	data <- sweep(data, 2, cm)
	data <- sweep(data, 2, sds, "/")
	
	#eliminate divide by zero
	wh0 <- which(sds == 0)
	if(length(wh0) > 0) {
        if(print == T) {
	        cat("Monitor ", i, ": No variability in ", colnames(data)[wh0], "\n")
        }
            
        data[, wh0] <- NA
	    data <- data[, -wh0]
	}
	
	list(data = data, wh0 = wh0)
	
}




getbstar <- function(stddata, nsources, dates) {
		
	#perform PCA
	pc1 <- prcomp(stddata)
	#get eigenvectors scaled by sdev
	rots <- sweep(pc1$rot, 2, pc1$sdev, "*")
	
	#apply varimax rotation
	vmax <- varimax(rots[, 1 : nsources], normalize = T)
	bstar1 <- as.matrix(vmax[[1]][1 : ncol(stddata), ])
	
	#inverse correlation matrix
	scor <- chol2inv(chol(cor(stddata)))
	bstar <- scor %*% bstar1
	scores <- as.matrix(stddata) %*% bstar
	
	rownames(scores) <- dates

	res <- list(scores = scores, bstar = bstar, vmax = vmax)
	res
}


	




#Function to get scores and profiles
getscoresprofs <- function(data, bstar, scores, tots, dates, mons = NULL,
    type = "apca") {

	scores1 <- scores
    
    if(tolower(type) == "mapca" & is.null(mons)) {
        stop("Need monitor list for mAPCA")
    }
	
	cc <- complete.cases(t(scores))
	scores <- scores[, cc]
	nsources <- ncol(scores)
	
	#get pollution-free day
	pc0 <- abzero(data, bstar[, cc])

	#compute absolute principal component scores
	apcs <- sweep(scores, 2, pc0)
	
	
	#get PM mass (either given, or rowsums)
	if (is.null(tots)) {
		tots <- rowSums(data)
	}
    
    
    #get scores
	temp <- regress1(apcs, tots, mons, type)
    leftover <- temp$xi0
    betas <- temp$coefs
	conc1 <- sweep(apcs, 2, betas, "*")
    
    #get profiles
	profs <- getprofs(data, conc1, nsources, mons, type)

	
	conc <- matrix(nrow= nrow(conc1), ncol = ncol(scores1))
	conc[, cc] <- conc1
	colnames(conc) <- paste0("source", seq(1, nsources))
    
    if(tolower(type) == "mapca") {
        conc <- data.frame(conc, mons)
        rn <- paste0(dates, mons)
    }else{
        rn <- dates
        
    }
	rownames(conc) <- rn
	
	
	res <- list(conc = conc, profs = profs, leftover = leftover, betas = betas)
	res
	}
	

regress1 <- function(x, y, mons = NULL, type = "apca") {
    
    
    #get source concentrations
    if(tolower(type) == "apca") {
        reg1 <- lm(y ~ x)
        #unexplained sources
        leftover <- reg1$coef[1]
        betas <- reg1$coef[-1]
  
    }else if(tolower(type) == "mapca") {
        reg1 <- try(lme(y ~ x, random = ~1 | mons), silent = T)
        
        if(class(reg1) == "try-error") {
            reg1 <- try(lme(y ~ x, random = ~1 | mons, 
                            control = lmeControl(opt = "optim")))
        }
        leftover <- reg1$coef[[1]][1]
        betas <- reg1$coef[[1]][-1]
        
    }
    
    list(xi0 = leftover, coefs = betas)
    
}
	
#Function to get profiles for APCA
getprofs <- function(data, conc1, nsources, mons, type) {

    #get source profiles

    profs <- matrix(nrow = nsources, ncol = ncol(data))
    for (i in 1 : ncol(data)) {
          
        lm1 <- try(regress1(conc1, data[, i], mons, type))
        
        if (class(lm1) != "try-error") {
            profs[, i] <- lm1$coefs
        }
            
    }	
    
    rownames(profs) <- paste0("source", seq(1, nsources))
    profs
}

#Get rotated pollution-free day
#data is uncentered data
#bstar is factor matrix
abzero <- function(data, bstar){
	
	#find column means and sd
	colm <- apply(data, 2, mean, na.rm = T)
	colsd <- apply(data, 2, sd, na.rm = T)
	
	#zero pollution day
	z0 <- -colm / colsd
	
	whNAN <- which(is.nan(z0))
	if(length(whNAN) > 0) {
		z0 <- z0[-whNAN]
	}
	
	#rotated into factor space
	p0 <- t(bstar) %*% z0

	p0	
}



#' @export
print.apca <- function(x) {
	cat("Call:\n")
	print(x$call)
	cat("Profiles:\n")
	print(x$prof)
}


#' @export
summary.apca <- function(x, ...) {
	
	ns <- ncol(x$conc)
	means <- apply(x$conc, 2, mean, ...)
	sd <- apply(x$conc, 2, sd, ...)
	meanssd <- rbind(means, sd)
	
	res <- list(nsources = ns, meanssd = meanssd,
		profs = round(x$prof, 2), call = match.call())
	class(res) <- "summary.apca"
	res
}

#' @export
print.summary.apca <- function(x) {
	cat("Call:\n")
	print(x$call)
	cat("\n")
	cat("Number of sources:", x$nsources, "\n")
	cat("Summary of source concentrations:\n")
	print(x$meanssd)
	cat("Source profiles:\n")
	print(x$profs)
}

#' @export
plot.apca <- function(x, plot = "prof", names = NULL, 
	dates1 = NULL) {

	if(is.null(names)) {
		names <- colnames(x$conc)
	}
	
	
	if(plot == "prof") {
		names 
		barplot(x$prof, beside = T, 
			legend.text = names)
	}else if(plot == "conc") {
		conc <- x$conc
		dates <- x$dates
		
		#select sources to plot
		if(length(names) > 9) {
			print("Randomly selected 9 sources")
			names <- sample(names, 9)
		}
		
		#either print subset (if # dates > 200)
		if(nrow(conc) > 200 & is.null(dates1)) {
			print("Plotting first 200 days")
			conc <- conc[1 : 200, ]
			dates <- dates[1: 200]
		#or dates are specified in function call
		}else if(!is.null(dates1)){
			whD <- which(dates %in% dates1)
			dates <- dates[whD]
			conc <- conc[whD, ]
			
		}
		
		
		par(mfrow = c(3, 3))
		for(i in 1 : length(names)) {
			plot(dates, conc[, names[i]], type = "l", axes = F, 
				xlab = "", ylab = "Concentration", main = names[i])
			axis.Date(1, dates, las = 2, format = "%m/%Y")
			axis(2)
			box()
			abline(h = 0, col = "blue", lty = 2)
			}
		
	}
	
}


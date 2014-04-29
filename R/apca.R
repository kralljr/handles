######
# Functions to perform APCA


apca <- function(x, ...) UseMethod("apca")


#function to perform apca
# data is chemical constituent data (days X # constituents)
#tots (optional) is vector of PM for each day
#nfactors is number of sources
apca.default <- function(data, tots = NULL,
	nsources = NULL, bstar1 = NULL, adjust = NULL,  
	mdl = NULL, ...){
		
	if(!is.null(adjust)) {
		adj1 <- adjust(data = data, mdl = mdl, 
			method = adjust, ...)
		dat <- adj1$dat
	}	
		
	if(class(dat[, 1]) != "Date") {
		stop("First column must be 'date'")
	}
		
	dates <- data[, 1]
	data <- data[, -1]	
		
	#standardize data (mean zero, sd 1)
	stddata <- stdize1(data, ...)
	
	#Use PC method to get number of sources
	if(is.null(nsources)) {
		pr1 <- prcomp(stddata, ...)
		nsources <- length(which(pr1$sdev > 1))
	}
		
	#get scores and rotation matrix bstar	
	temp <- getbstar(stddata = stddata, nsources = nsources, 
		dates = dates, ...)
	scores <- temp$scores
	bstar <- temp$bstar
	vmax <- temp$vmax

	#use regressions to obtain source profiles and source concentrations
	apca <- getscoresprofs(data = data, bstar = bstar, 
		scores = scores, tots = tots, dates = dates, ...) 

	apca$vmax <- vmax
	apca$dates <- dates
	apca$nsources <- nsources
	apca$call <- match.call()
	
	apca$adjust <- adj1

	class(apca) <- "apca"
	apca
}





	
	
#### function to standardize data: mean 0, var 1
# dat is data matrix T (ndays) X P (nconstituents)
stdize1 <- function(data, ...) {
	
	#get means and sd
	cm <- colMeans(data, na.rm = T, ...)
	sds <- apply(data, 2, sd, na.rm = T, ...)
	
	#standardize
	data <- sweep(data, 2, cm, ...)
	data <- sweep(data, 2, sds, "/", ...)
	
	#eliminate divide by zero
	wh0 <- which(sds == 0)
	if(length(wh0) > 0) {
		data <- data[, -wh0]
		}

	data
	
}




#get factor rotation and scores
#datasc is centered and scaled data
#nsource is number of sources
getbstar <- function(stddata, nsources, dates, ...) {
		
	#perform PCA
	pc1 <- prcomp(stddata, ...)
	#get eigenvectors scaled by sdev
	rots <- sweep(pc1$rot, 2, pc1$sdev, "*", ...)
	
	#apply varimax rotation
	vmax <- varimax(rots[, 1 : nsources], normalize = T, ...)
	bstar1 <- as.matrix(vmax[[1]][1 : ncol(stddata), ], ...)
	
	#get uncorrelated factors and scores

	#inverse correlation matrix
	scor <- chol2inv(chol(cor(stddata, ...), ...), ...)
	bstar <- scor %*% bstar1
	scores <- as.matrix(stddata, ...) %*% bstar
	
	rownames(scores) <- dates

	res <- list(scores = scores, bstar = bstar, vmax = vmax)
	res
}


	





getscoresprofs <- function(data, bstar, scores, tots, dates, ...) {

	scores1 <- scores
	
	cc <- complete.cases(t(scores))
	scores <- scores[, cc]
	nsources <- ncol(scores)
	
	#get pollution-free day
	pc0 <- abzero(data, bstar[, cc], ...)

	#compute absolute principal component scores
	apcs <- sweep(scores, 2, pc0)
	
	
	#get PM mass (either given, or rowsums)
	if (is.null(tots)) {
		tots <- rowSums(data)
	}
		
	#get source concentrations	
	reg1 <- lm(tots ~ apcs, ...)
	
	
	#unexplained sources
	leftover <- reg1$coef[1]
	betas <- reg1$coef[-1]
	conc1 <- sweep(apcs, 2, betas, "*", ...)


	#get source profiles
	
	
	Xmat <- cbind(rep(1, nrow(conc1)), conc1)
	#get (X'X)^(-1)
	s1 <- chol(t(Xmat) %*% Xmat)
	s1 <- try(chol2inv(Xmat))
	
	if (class(s1) != "try-error"){
		
		betas <- s1 %*% t(Xmat) %*% as.matrix(data)
		profs <- betas[-1, ]  
		
	} else {
		
		profs <- matrix(nrow = nsources, ncol = ncol(data))
		for (i in 1 : ncol(data)) {
			lm1 <- try(lm(data[, i] ~ conc1))
			
			if (class(lm1) != "try-error") {
				profs[, i] <- lm1$coef[-1]
			}

		}
		
	}	
	
	conc <- matrix(nrow= nrow(conc1), ncol = ncol(scores1))
	conc[, cc] <- conc1
	colnames(conc) <- paste0("source", seq(1, nsources))
	rownames(conc) <- dates
	rownames(profs) <- paste0("source", seq(1, nsources))
	
	
	res <- list(conc = conc, profs = profs, leftover = leftover)
	res
	}
	
	
	


#Get rotated pollution-free day
#data is uncentered data
#bstar is factor matrix
abzero <- function(data, bstar, ...){
	
	#find column means and sd
	colm <- apply(data, 2, mean, na.rm = T, ...)
	colsd <- apply(data, 2, sd, na.rm = T, ...)
	
	#zero pollution day
	z0 <- -colm / colsd
	
	whNAN <- which(is.nan(z0), ...)
	if(length(whNAN) > 0) {
		z0 <- z0[-whNAN]
	}
	
	#rotated into factor space
	p0 <- t(bstar) %*% z0

	p0	
}


	


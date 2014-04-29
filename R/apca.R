######
# Functions to perform APCA





#function to perform apca
# data is chemical constituent data (days X # constituents)
#tots (optional) is vector of PM for each day
#nfactors is number of sources
apca <- function(data, tots = NULL,
	nfactors=FALSE, bstar1 = NULL){
		
	#get mean zero, sd1 data
	datasc <- stdize1(data)[[1]]
	
	#if nfactors == false, use PC method to get # sources
	if(nfactors==FALSE){

		pcgetnc <- prcomp(datasc)
		nsource <- length(which(pcgetnc$sdev > 1))
		
	}else{
		nsource <- nfactors
		}
		
		
	#get scores and rotation matrix bstar	
	if(is.null(bstar1)) {
		temp <- getbstar(datasc, nsource)
		pc <- temp[[1]]
		bstar <- temp[[2]]
		vmax <- temp[[3]]
		vls <- temp[[4]]
	}else{
		scor <- chol2inv(chol(cor(datasc)))
		bstar <- scor %*% bstar1
		pc <- as.matrix(datasc) %*% bstar
		rots <- NULL
		}

	#use regressions to obtain source profiles and source concentrations
	tout <- getscoresprofs(data, bstar, pc, tots) 


	scores <- tout[[1]]
	profs <- tout[[2]]
	xi0 <- tout[[3]]
	betas <- tout[[5]]

	
	list(scores, bstar, xi0, t(profs), vmax, betas, vls)
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


#get factor rotation and scores
#datasc is centered and scaled data
#nsource is number of sources
getbstar <- function(datasc, nsource) {
		
		#perform PCA
		pc1 <- prcomp(datasc)
		#get eigenvectors scaled by sdev
		rots <- sweep(pc1$rot, 2, pc1$sdev, "*")
		
		#apply varimax rotation
		vmax1 <- varimax(rots[, 1 : nsource], normalize = T)
		bstar1 <- as.matrix(vmax1[[1]][1 : ncol(datasc), ])
		
		#get uncorrelated factors and scores
		# bstar <- solve(cor(datasc)) 
		# browser()
		scor <- chol2inv(chol(cor(datasc)))
		bstar <- scor %*% bstar1
		scores <- as.matrix(datasc) %*% bstar


	list(scores, bstar, vmax1, bstar1)	
}


	
	
	
#### function to standardize data: mean 0, var 1
# dat is data matrix T (ndays) X P (nconstituents)
stdize1 <- function(dattemp) {
	
	#get means and sd
	cm <- colMeans(dattemp, na.rm = T)
	sds <- apply(dattemp, 2, sd, na.rm = T)
	
	#standardize
	dattemp <- sweep(dattemp, 2, cm)
	dattemp <- sweep(dattemp, 2, sds, "/")
	
	#eliminate divide by zero
	wh0 <- which(sds == 0)
	if(length(wh0) > 0) {
		# print(c("wh0", wh0))
		dattemp[, wh0] <- NA
		dattemp <- dattemp[, -wh0]
		}

	list(dattemp, wh0)
	
}




getscoresprofs <- function(data, bstar, pc, tots) {
	
	

	pcold <- pc
	
	cc <- complete.cases(t(pc))
	pc <- pc[, cc]
	nsource <- ncol(pc)
	#get pollution-free day
	pc0 <- abzero(data, bstar[, cc])

	#compute absolute principal component scores
	pc0 <- matrix(rep(pc0, each=nrow(pc)), byrow = FALSE,
		ncol = nsource)
	apca <- pc - pc0
	
	
	#get PM mass (either given, or rowsums)
	if (is.null(tots)) {
		totmass <- rowSums(data)
	} else {
		totmass <- tots
		}
		
	#get source concentrations	
	reg1 <- lm(totmass ~ apca)
	#unexplained sources
	xi0 <- reg1$coef[1]
	betas <- reg1$coef[-1]
	scores <- matrix(rep(betas, each = nrow(data)),
		byrow = FALSE, ncol = nsource) * apca 
			
	#get source profiles
	scores2 <- cbind(rep(1, nrow(scores)), scores)
	s1 <- chol(t(scores2) %*% scores2)
	s1 <- try(chol2inv(s1)) # t(ncons+1 X nsource)
	prof <- matrix(ncol = nsource, nrow = ncol(data))
	if (class(s1) != "try-error" | length(which(is.na(s1))) > 0){
		betas <- s1 %*% t(scores2) %*% as.matrix(data)
		profs <- betas[-1, ]  #nsource X ncons
		
	} else {
		profs <- matrix(nrow = nsource, ncol = ncol(data))
		for (i in 1 : ncol(data)) {
			
			lm1 <- try(lm(data[, i] ~ scores))
			
			if (class(lm1) != "try-error") {
				profs[, i] <- lm1$coef[-1]
			}
			
			
		}
	}	
	
	scores2 <- matrix(nrow= nrow(scores), ncol = ncol(pcold))
	scores2[, cc] <- scores
	
	list(scores2, profs, xi0, cc, betas)
	}


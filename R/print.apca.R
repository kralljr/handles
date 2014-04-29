####
# Functions for printing APCA results


print.apca <- function(x, ...) {
	cat("Call:\n")
	print(x$call)
	cat("Profiles:\n")
	print(x$prof)
}

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


print.summary.apca <- function(x, ...) {
	cat("Call:\n")
	print(x$call)
	cat("\n")
	cat("Number of sources:", x$nsources, "\n")
	cat("Summary of source concentrations:\n")
	print(x$meanssd)
	cat("Source profiles:\n")
	print(x$profs)
}


plot.apca <- function(x, plot = "prof", names = NULL, 
	dates1 = NULL, ...) {

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
		if(nrow(conc) > 200 & is.null(dates)) {
			print("Plotting first 200 days")
			conc <- conc[1 : 200, ]
			dates <- dates[1: 200]
		#or dates are specified in function call
		}else{
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
			}
		
	}
	
}
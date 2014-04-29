####
# Functions for printing APCA results


print.apca <- function(x, ...) {
	print("Call:\n")
	print(x$call)
	cat("Profiles:\n")
	print(x$prof)
}

summary.apca <- function(x, ...) {
	cat("Number of sources:", x$nf, "\n")
	cat("Summary of source concentrations:\n")
	print(head(x$conc))
}


plot.apca <- function(x, ...) {
	
	
}
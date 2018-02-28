#' exac
#' 
#' Query ExAC database for Consequence, Allele Frequency, PolyPhen Score and
#' gene SYMBOL
#'
#' @param var A single variant in the format chr-pos-ref-alt
#' 
#' @return A vector of the four items to look up 
#' @author Alexander Hopkins
#' @export

#a function to call the shell script which does the ExAC lookup
exac <- function(var){
	out <- system(paste('sh exac_lookup.sh -cfps',var), intern=T)
	#clean up the output (NA instead of empty or 'null')
	out[out==''] <- NA
	out[out=='null'] <- NA
	#return
	out
}


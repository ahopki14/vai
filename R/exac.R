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
	#locate the lookup script in the package and call it
	if(!length(out)>0){warning('Nothing generated from the query,
				  make sure that curl is installed and 
				  the exac_lookup.sh script is executable'
				  )
	}
	#return
	out
}


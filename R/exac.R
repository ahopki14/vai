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
	path <- system.file('exec','exac_lookup.sh',package='vai')
	out <- system(paste('sh',path,'-cfps',var), intern=T)
	#clean up the output (NA instead of empty or 'null')
	out[out==''] <- NA
	out[out=='null'] <- NA
	#check that something was created
	if(!length(out)>0){warning('Nothing generated from the query,
				  make sure that curl is installed and 
				  the exac_lookup.sh script is executable'
				  )
	}
	#return
	out
}

#make some change


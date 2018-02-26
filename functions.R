
#a function to call the shell script which does the ExAC lookup
exac <- function(var){
	out <- system(paste('sh exac_lookup.sh -cfp',var), intern=T)
	#clean up the output (NA instead of empty or 'null')
	out[out==''] <- NA
	out[out=='null'] <- NA
	#return
	out
}



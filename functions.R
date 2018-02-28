#a function to call the shell script which does the ExAC lookup
exac <- function(var){
	out <- system(paste('sh exac_lookup.sh -cfps',var), intern=T)
	#clean up the output (NA instead of empty or 'null')
	out[out==''] <- NA
	out[out=='null'] <- NA
	#return
	out
}

#To solve the 'minimal representation problem' I borrowed an idea from 
#http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
# which was linked from the ExAC documentation

simplify_variant <- function(x){
	stopifnot(class(x)=='ExpandedVCF' & length(x)==1)
	nm <- rownames(x)
	#break up into CHR, POS REF and ALT
	sp <- strsplit(nm,':|_')[[1]]
	sp <- c(sp[-3],strsplit(sp[3],'/')[[1]])
	POS <- as.numeric(sp[2])
	REF <- strsplit(as.character(ref(x)),'')[[1]]
	ALT <- strsplit(as.character(alt(x)),'')[[1]]
	#trim off any matching bases at the end
	while(REF[length(REF)]==ALT[length(ALT)] & min(length(REF),length(ALT))>1){
		REF <- REF[-length(REF)]
		ALT <- ALT[-length(ALT)]
	}
	#then trim any at the beginning, adjusting the position up accordingly
	while(REF[1]==ALT[1] & min(length(REF),length(ALT))>1){
		REF <- REF[-1]
		ALT <- ALT[-1]
		POS <- POS+1
	}
	#fix the record
	REF <- paste0(REF, collapse='')
	ALT <- paste0(ALT, collapse='')
	out <- x
	ref(out) <- DNAStringSet(REF)
	alt(out) <- DNAStringSet(ALT)
	#rownames(x) <- paste0(sp[1],':',POS,'_',REF,'/',ALT)
	names(out@rowRanges)[1] <- paste0(sp[1],':',POS,'_',REF,'/',ALT)
	#return x
	out
}


#a function to call the shell script which does the ExAC lookup
exac <- function(var){
	out <- system(paste('sh exac_lookup.sh -cfps',var), intern=T)
	#clean up the output (NA instead of empty or 'null')
	out[out==''] <- NA
	out[out=='null'] <- NA
	#return
	out
}

#a function to decompose multi-allelic variants into individual lines
#this is no longer needed when using expand(vcf)
decomp_vcf <- function(vcf){
	stopifnot(class(vcf)=='CollapsedVCF')
	out <- vcf[0] #makes an empty VCF with the right INFO, etc. 
	for(a in nrow(vcf)){
		record <- vcf[a]
		tmp <- record[c(1,1)] #duplicate the record
		geno(tmp)$GT <- matrix(c('0/1','1/0')) 
		#this is ugly, but they need to be formatted as a list
		geno(tmp)$AO[1] <- list(unlist(geno(tmp)$AO[1])[1])
		geno(tmp)$AO[2] <- list(unlist(geno(tmp)$AO[2])[2])
		info(tmp)$AF[1] <- NumericList(info(tmp)$AF[[1]][1])
		info(tmp)$AF[2] <- NumericList(info(tmp)$AF[[2]][2])
		fixed(tmp)$ALT[1] <- DNAStringSetList(fixed(tmp)$ALT[[1]][1])
		fixed(tmp)$ALT[2] <- DNAStringSetList(fixed(tmp)$ALT[[2]][2])
		#Simplify the variants
		ALT <- as.character(fixed(tmp)$ALT[[2]])
		REF <- as.character(fixed(tmp)$REF[[2]])
		fixed <- simplify_variant(rownames(record)[1],REF,ALT)
		#overwrite the existing fields with the simplified variant
		fixed(tmp)$REF[2] <- DNAStringSet(fixed$REF)
		fixed(tmp)$ALT[2] <- DNAStringSetList(fixed$ALT)
		#fix the rownames
		rownames(tmp) <- c(rownames(record)[1], fixed$name)
		#return the VCF object with 2 rows
		out <- rbind(out,tmp)
	}
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


#a function to call the shell script which does the ExAC lookup
exac <- function(var){
	out <- system(paste('sh exac_lookup.sh -cfp',var), intern=T)
	#clean up the output (NA instead of empty or 'null')
	out[out==''] <- NA
	out[out=='null'] <- NA
	#return
	out
}

#a function to decompose multi-allelic variants into individual lines
decomp_vcf <- function(x){
	stopifnot(class(x)=='CollapsedVCF' & length(x)==1)
	tmp <- x[c(1,1)] #duplicate the record
	geno(tmp)$GT <- matrix(c('0/1','1/0')) 
	#this is ugly, but they need to be formatted as a list
	geno(tmp)$AO[1] <- list(unlist(geno(tmp)$AO[1])[1])
	geno(tmp)$AO[2] <- list(unlist(geno(tmp)$AO[2])[2])
	info(tmp)$AF[1] <- NumericList(info(tmp)$AF[[1]][1])
	info(tmp)$AF[2] <- NumericList(info(tmp)$AF[[2]][2])
	fixed(tmp)$ALT[1] <- DNAStringSetList(fixed(tmp)$ALT[[1]][1])
	fixed(tmp)$ALT[2] <- DNAStringSetList(fixed(tmp)$ALT[[2]][2])
	#last, fix the rownames, needed later for ExAC lookup
	rownames(info(tmp)) <- c('1','2') #needed to establish rownames
	rownames(info(tmp))[1] <- rownames(x)
	rownames(info(tmp))[2] <- gsub(pattern='\\/[ACTG]*',
				       replacement=paste0('/',
						      as.character(fixed(tmp)$ALT[[2]])
						      ),
				       x=rownames(x)
				       )
}

#the problem with the decomp_vcf function is that the variant names are not
#sufficiently simplified
simplify_name <- function(nm){
	#break up into CHR, POS REF and ALT
	sp <- strsplit(nm,':|_')[[1]]
	sp <- c(sp[-3],strsplit(sp[3],'/')[[1]])
	POS <- as.numeric(sp[2])
	#Identify which positions are mutated
	REF <- strsplit(sp[3],'')[[1]]
	ALT <- strsplit(sp[4],'')[[1]]
	mutated <- !(REF==ALT)
	#trim off any matching bases at the end
	while(REF[length(REF)]==ALT[length(ALT)]){
		REF <- REF[-length(REF)]
		ALT <- ALT[-length(ALT)]
	}
	#then trim any at the beginning, adjusting the position up accordingly
	while(REF[1]==ALT[1]){
		REF <- REF[-1]
		ALT <- ALT[-1]
		POS <- POS+1
	}
	#return the corrected name
	paste0(sp[1],':',POS,'_',REF,'/',ALT)
}

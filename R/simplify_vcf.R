#a function to apply the simplify_variant function to all
#necessary records of a vcf
simplify_vcf <- function(vcf){
	stopifnot(class(vcf)=='ExpandedVCF')
	#identify rows which require fixes and fix them
	#  the variant names, reference and alternate alleles and 
	#  position are all fixed 
	w <- grep('2', geno(vcf)$GT)
	fixed_records  <- sapply(vcf[w], simplify_variant)

	#sanity check that all records were checked
	stopifnot(length(fixed_records)==length(w))

	#put the new records back in the original vcf data file
	for(a in seq_along(fixed_records)){
		vcf[w[a]] <- fixed_records[[a]]
	}
	vcf
}

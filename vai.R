#-------------Load libraries and data --------------------------------------
library(VariantAnnotation)
library(dplyr)

v_total <- readVcf(file='./data/coding_challenge_final.vcf')
source('functions.R')

#------------- Fix multi-allelic variants ----------------------------------

#expand VCF to one allele per line 
v_expanded <- expand(v_total, row.names=T)

#identify rows which require fixes and fix them
#  the variant names, reference and alternate alleles and 
#  position are all fixed 
w <- grep('2', geno(v_expanded)$GT)
fixed_records  <- sapply(v_expanded[w], simplify_variant)

#sanity check that all records were checked
stopifnot(length(fixed_records)==length(w))

#put the new records back in the original vcf data file
for(a in seq_along(fixed_records)){
	v_expanded[w[a]] <- fixed_records[[a]]
}

#----------- Extract info from VCF ---------------------------------------

#get necessary columns from the INFO
info_data <- info(v_expanded)[,c('TYPE','AF','DP')]


#get necessary columns from the SAMPLE field
sample_data <- DataFrame(GT=geno(v_expanded)$GT,
			  AO=unlist(geno(v_expanded)$AO),
			  RO=geno(v_expanded)$RO
			  )
colnames(sample_data) <- c('GT','AO','RO') #not sure why this doesn't default

#sanity check that the two DataFrames are in the same order
stopifnot(all(rownames(sample_data)==rownames(info_data)))

#put them together
data <- cbind(info_data, sample_data)

#make a column for percentage of reads supporting alternate
data$PctAlt <- round(data$AO/data$DP,2)


#----------- Query ExAC for other info ---------------------------------------

#Format the variant names for ExAC
nms <- rownames(data)
nms <- gsub('chr','',nms)
nms <- gsub(':|_|\\/','-',nms)
data$VariantNames <- nms

# Get the ExAC info for all 
start <- proc.time()
	out <- mapply(data$VariantNames, FUN=exac)#this ran in ~8 minutes
	#check the time to query
	t <- round((proc.time()-start)[['elapsed']]/60,2)
cat(dim(out)[2],'queries in',t,'minutes\n') 

#Clean up the ExAC data
out <- DataFrame(t(out))
names(out) <- c('Consequence', 'PopAlleleFreq','PolyPhen','GeneName')
out$Consequence <- as.character(out$Consequence)
out$GeneName <- as.character(out$GeneName)
out$PolyPhen <- as.numeric(as.character(out$PolyPhen))
out$PopAlleleFreq <- as.numeric(as.character(out$PopAlleleFreq))
out$PopAlleleFreq[is.na(out$PopAlleleFreq)] <- 0
out$Consequence[is.na(out$Consequence)] <- 'Unknown'

#sanity check the order
stopifnot(all(rownames(out)==data$nms))


#----------- Combine VCF and ExAC Data and save----------------------------
data <- cbind(data,DataFrame(out))

#make names more readable for an Investigator
data <- rename(as.data.frame(data),
	       Type=TYPE,
	       AllelFreq=AF,
	       ReadDepth=DP,
	       Genotype=GT,
	       AltReads=AO,
	       RefReads=RO
	       )

#order by type and PolyPhen score and save
o <- order(data$Type,data$PolyPhen,decreasing=T)
write.csv(data[o,],file='table.csv', row.names=F)


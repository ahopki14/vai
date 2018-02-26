library(VariantAnnotation)

v_total <- readVcf(file='./data/coding_challenge_final.vcf')

#decompose multi-allelic variants 
v_expanded <- expand(v_total, row.names=T)

#identify rows which require fixed rownames and fix them
w <- grep('2', geno(v_expanded)$GT)
for(a in w){
	v_expanded[a] <- simplify_variant(v_expanded[a])
}


#get necessary columns from the INFO
info_data <- info(v_expanded)[,c('TYPE','DPRA','AF','DP')]


#get necessary columns from the SAMPLE field
sample_data <- DataFrame(GT=geno(v_singles)$GT,
			  AO=unlist(geno(v_singles)$AO),
			  RO=geno(v_singles)$RO
			  )
colnames(sample_data) <- c('GT','AO','RO') #not sure why this doesn't default

#sanity check that the two DataFrames are in the same order
stopifnot(all(rownames(sample_data)==rownames(info_data)))

#put them together
data <- cbind(info_data, sample_data)


#Format the variant names for ExAC
nms <- rownames(data)
nms <- gsub('chr','',nms)
nms <- gsub(':|_|\\/','-',nms)
data$nms <- nms

# Get the ExAC info for all 
start <- proc.time()
out <- mapply(data$nms, FUN=exac)
#check the time to query
t <- round((proc.time()-start)[['elapsed']]/60,2)
cat(dim(out)[2],'queries in',t,'minutes\n')

out <- as.data.frame(t(out))
names(out) <- c('Consequence', 'Pop_Allele_Freq','PolyPhen')
out$PolyPhen <- as.numeric(out$PolyPhen)
out$Pop_Allele_Freq <- as.numeric(out$Pop_Allele_Freq)

#sanity check the order
stopifnot(all(rownames(out)==data$nms))

data <- cbind(data,out)
write.csv(data,file='out.csv')


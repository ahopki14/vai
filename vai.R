library(VariantAnnotation)

v_total <- readVcf(file='./data/coding_challenge_final.vcf')

#for now, set aside multiallelic variants
w <- grep('2',geno(v_total)$GT)
v_singles <- v_total[-w,]

#get necessary columns from the INFO
info_data <- info(v_singles)[,c('TYPE','DPRA','AF','DP')]


#get necessary columns from the SAMPLE field
sample_data <- DataFrame(GT=geno(v_singles)$GT,
			  AO=unlist(geno(v_singles)$AO),
			  RO=geno(v_singles)$RO
			  )
colnames(sample_data) <- c('GT','AO','RO') #not sure why this doesn't default

#sanity check that the two DataFrames are in the same order
stopifnot(all(rownames(sample_data)==rownames(i)))

#put them together
data <- cbind(info_data, sample_data)


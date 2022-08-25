library(dplyr)
library(tidyr)
library(survival)
library(parallel)
library(data.table)
library(RhpcBLASctl)



blas_set_num_threads(5)

infile= snakemake@input[[1]]
phenofile= snakemake@input[[2]]
ID= 'SentrixID_1' #'MOR_PID' # ID name
outfile= snakemake@output[[1]]
time_t= 'SVLEN_UL_DG' # time variable name
outcome= 'spont' # event variable name
covars= c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PARITY0')

options(stringsAsFactors=FALSE)

pheno= fread(phenofile)


pheno= pheno[!duplicated(pheno$IID), ]


dataChunk= fread(infile, header= T)
if (nrow(dataChunk)==0){
cat(NULL, file= snakemake@output[[1]])
cat(0, file= snakemake@output[[2]])
} else {
	dataChunk= dataChunk[rowSums(dataChunk[, 3:ncol(dataChunk)])>20, ]
        if (nrow(dataChunk)==0){
	cat(NULL, file= snakemake@output[[1]])
cat(0, file= snakemake@output[[2]])
} else {
	genvars= paste(dataChunk$CHR, dataChunk$segment, sep=':')
	dataChunk= subset( dataChunk, select = -c(CHR,segment))
	
        dataChunk= as.data.frame(t(dataChunk))
	
        dataChunk$id= gsub('X','',rownames(dataChunk))
	names(dataChunk)[1:length(genvars)]= genvars
        geno= inner_join(pheno, dataChunk, by= c('IID' = 'id'))
	df_list= lapply(names(geno[,-c(1:dim(pheno)[2])]), function(snp){cox_coef= survreg(Surv( geno$SVLEN_UL_DG, geno$spont)~ geno[,snp] + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6 + geno$PC7 + geno$PC8 + geno$PC9 + geno$PC10 + geno$cohort + geno$PARITY0 + geno$FKB, na.action = na.omit, dist= 'weibull')
	n= summary(cox_coef)$n
       coef = summary(cox_coef)$table[2, 1]
	freq= mean(geno[, snp], na.rm=T)
	loglikF= cox_coef$loglik[2]
       sd= summary(cox_coef)$table[2,2]
       pvalue= summary(cox_coef)$table[2, 4]
#	txt = sprintf( "%s\t%e\t%e\t%e\t%e\t%e\n", snp, n, coef, sd, pvalue, loglikF)
#cat(txt, file= outfile, append= T)
	return(list(snp, n, freq, coef, sd, pvalue, loglikF))

}
)

z= do.call('rbind', df_list)
write.table(z, outfile, sep= '\t', row.names=F, col.names=F, quote= F)

#geno= as.matrix(geno[,-c(1:dim(pheno)[2])])
geno= geno[,-c(1:dim(pheno)[2])]
geno= geno[,!apply(geno, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
pca.fit= prcomp(geno, scale=TRUE, center= TRUE)

eff= sum(summary(pca.fit)$importance[3,]<0.995)

cat(eff, file= snakemake@output[[2]],sep= '\n')
}
}

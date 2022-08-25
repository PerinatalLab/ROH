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

options(stringsAsFactors=FALSE)

pheno= fread(phenofile)


pheno= pheno[!duplicated(pheno$IID), ]


dataChunk= fread(infile, sep= "\t")
        genes= dataChunk$gene
        dataChunk= subset(dataChunk, select = -c(gene))


        dataChunk= as.data.frame(t(dataChunk))
        dataChunk$id= gsub('X','',rownames(dataChunk))
        names(dataChunk)[1:length(genes)]= genes
        geno= inner_join(pheno, dataChunk, by= c('IID' = 'id'))
        cox_coef= lapply(names(geno[,-c(1:dim(pheno)[2]), drop=F]), function(snp){cox_coef= survreg(Surv( geno$SVLEN_UL_DG, geno$spont)~ geno[, snp] + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6 + geno$PC7 + geno$PC8 + geno$PC9 + geno$PC10 + geno$cohort + geno$PARITY0 + geno$FKB, na.action = na.omit, dist= 'weibull')
       n= summary(cox_coef)$n
       coef = summary( cox_coef)$coefficients[2]
	loglikF= cox_coef$loglik[2]
       sd= summary(cox_coef)$table[2,2]
	freq= mean(geno[, snp], na.rm=T)
       pvalue= summary(cox_coef)$table[2, 4]
        txt = sprintf( "%s\t%e\t%e\t%e\t%e\t%e\t%e\n", snp, n, freq, coef, sd, pvalue, loglikF)
cat(txt, file= outfile, append= T)
}
)

geno= geno[,-c(1:dim(pheno)[2])]

if (is.vector(geno)){
eff= 1
} else {
geno= geno[,!apply(geno, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
pca.fit= tryCatch(prcomp(geno, scale=TRUE, center= TRUE), warning = function(w) {print('Warning')}, error = function(e) { prcomp(geno, scale=FALSE, center= TRUE)})

eff= sum(summary(pca.fit)$importance[3,]<0.995)
}

cat(eff, file= snakemake@output[[2]],sep= '\n')



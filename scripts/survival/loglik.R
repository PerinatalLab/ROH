library(dplyr)
library(tidyr)
library(survival)
library(parallel)
library(data.table)
library(RhpcBLASctl)



blas_set_num_threads(2)

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

# read pheno file; each row is one participant, and column represents one variable
print(infile)

colnames= readLines(gzfile(infile), n= 1)
colnames= unlist(strsplit(colnames, '\t'))

pheno= filter(pheno, IID %in% colnames)

cox_coef= survreg(Surv(SVLEN_UL_DG, spont)~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +  PC7 + PC8 + PC9 + PC10 + cohort + PARITY0, data= pheno, na.action = na.omit, dist= 'weibull')
loglikF= cox_coef$loglik[2]
txt = sprintf( "%e\t\n", loglikF)
cat(txt, file= outfile, append= T)


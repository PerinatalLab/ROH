library(dplyr)
library(tidyr)
library(survival)
library(data.table)
library(RhpcBLASctl)

blas_set_num_threads(4)


infile= snakemake@input[[1]]

phenofile= snakemake@input[[2]]
ID= 'SentrixID_1' #'MOR_PID' # ID name
outfile= snakemake@output[[1]]
time_t= 'SVLEN_UL_DG' # time variable name
outcome= 'spont' # event variable name
covars= c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PARITY0')

options(stringsAsFactors=FALSE)

pheno= fread(phenofile)

if (file.size(infile)<=4) {
d= data.frame()
write(d, outfile)
quit(status=0)
} else {




pheno= pheno[!duplicated(pheno$IID), ]


ids= readLines(snakemake@input[[3]])

colnames= c(list('chr', 'pos', 'ref', 'eff'), ids)
colnames= unlist(colnames, use.names=FALSE)


funk= function(block.text){
	dataChunk= fread(text=block.text, sep="\t", col.names= colnames)
	if (ncol(dataChunk)==0) return(NULL)
	df= as.matrix(dataChunk[, 5:ncol(dataChunk)], ncol= length(5:ncol(dataChunk)), nrow= nrow(dataChunk))
	df[(df== '1|0') | (df== '1/0')]= 1
	df[(df== '0|1') | (df== '0/1')]= 1
	df[(df== '0|0') | (df== '0/0')]= 0
	df[(df== '1|1') | (df== '1/1')]= 2

	df= matrix(as.numeric(df), ncol= length(5:ncol(dataChunk)), nrow= nrow(dataChunk))
	eaf= rowSums(df) / (rowSums(!is.na(df)) * 2)
	
	df[which(eaf>0.5), ]= abs(df[which(eaf>0.5), ] -2)
	dataChunk= cbind(dataChunk[,1:4], df)
	dataChunk[which(eaf>0.5) + 4, c(3,4)]= dataChunk[which(eaf>0.5) + 4, c(4,3)]
	names(dataChunk)= colnames
	genvars= paste(dataChunk$chr, dataChunk$pos, dataChunk$ref, dataChunk$eff, sep=':')
	dataChunk= subset(dataChunk, select = -c(chr, pos, ref, eff))

	dataChunk= as.data.frame(t(dataChunk))
	names(dataChunk)[1:length(genvars)]= genvars
	dataChunk= dataChunk[,!apply(dataChunk, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
	if (!is.data.frame(dataChunk)) return(NULL)
        dataChunk$id= gsub('X', '', rownames(dataChunk)) 
        geno= inner_join(pheno, dataChunk, by= c('IID' = 'id'))

        cox_coef= lapply(names(geno[,-c(1:ncol(pheno))]), function(snp){
	hom= sum(geno[, snp]== 2)
	het= sum(geno[, snp]== 1)
	rec_m= geno[, snp]
	rec_m= ifelse(rec_m==2, 1, 0)
	dom_m= geno[, snp]
	dom_m= ifelse(dom_m>0, 1, 0)
	cox_coef= tryCatch(survreg(Surv(geno$SVLEN_UL_DG, geno$spont)~ rec_m + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6 + geno$PC7 + geno$PC8 + geno$PC9 + geno$PC10 + geno$cohort + geno$FKB, na.action = na.omit, dist= 'weibull'), warning = function(cond) {NA}, error = function(cond) {NA})
        if (is.na(cox_coef)) {
        ref= NA
	eff= NA
	n= NA
	coef=NA
	sd= NA
	pvalue= NA
#	correlation= NA 
#	corr_pvalue= NA
} else {
	ref= unlist(strsplit(snp, ':'))[[3]]
	eff= unlist(strsplit(snp, ':'))[[4]]
	snp= ifelse(unlist(strsplit(snp, ':'))[[3]] > unlist(strsplit(snp, ':'))[[4]], paste(unlist(strsplit(snp, ':'))[[1]], unlist(strsplit(snp, ':'))[[2]], unlist(strsplit(snp, ':'))[[4]], unlist(strsplit(snp, ':'))[[3]], sep=':'), snp)
        coef= summary( cox_coef)$coefficients[2]
        sd= summary(cox_coef)$table[2,2]
        n= summary(cox_coef)$n
        pvalue= summary(cox_coef)$table[2,4]
}
	txt = sprintf( "%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\n", snp, ref, eff, 'recessive', n, het, hom, coef, sd, pvalue)
cat(txt, file= outfile, append= T)

#Dominant
cox_coef= tryCatch(survreg(Surv(geno$SVLEN_UL_DG, geno$spont)~ dom_m + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6 + geno$PC7 + geno$PC8 + geno$PC9 + geno$PC10 + geno$cohort + geno$FKB, na.action = na.omit, dist= 'weibull'), warning = function(cond) {NA}, error = function(cond) {NA})
        if (is.na(cox_coef)) {
        ref= NA
        eff= NA
        n= NA
        coef=NA
        sd= NA
        pvalue= NA
} else {
        ref= unlist(strsplit(snp, ':'))[[3]]
        eff= unlist(strsplit(snp, ':'))[[4]]
        snp= ifelse(unlist(strsplit(snp, ':'))[[3]] > unlist(strsplit(snp, ':'))[[4]], paste(unlist(strsplit(snp, ':'))[[1]], unlist(strsplit(snp, ':'))[[2]], unlist(strsplit(snp, ':'))[[4]], unlist(strsplit(snp, ':'))[[3]], sep=':'), snp)
        coef= summary( cox_coef)$coefficients[2]
        sd= summary(cox_coef)$table[2,2]
        n= summary(cox_coef)$n
        pvalue= summary(cox_coef)$table[2,4]
}
        txt = sprintf( "%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\n", snp, ref, eff, 'dominant', n, het, hom, coef, sd, pvalue)
cat(txt, file= outfile, append= T)
# Additive
cox_coef= tryCatch(survreg(Surv(geno$SVLEN_UL_DG, geno$spont)~ geno[, snp] + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6 + geno$PC7 + geno$PC8 + geno$PC9 + geno$PC10 + geno$cohort + geno$FKB, na.action = na.omit, dist= 'weibull'), warning = function(cond) {NA}, error = function(cond) {NA})
        if (is.na(cox_coef)) {
        ref= NA
        eff= NA
        n= NA
        coef=NA
        sd= NA
        pvalue= NA
} else {
        ref= unlist(strsplit(snp, ':'))[[3]]
        eff= unlist(strsplit(snp, ':'))[[4]]
        snp= ifelse(unlist(strsplit(snp, ':'))[[3]] > unlist(strsplit(snp, ':'))[[4]], paste(unlist(strsplit(snp, ':'))[[1]], unlist(strsplit(snp, ':'))[[2]], unlist(strsplit(snp, ':'))[[4]], unlist(strsplit(snp, ':'))[[3]], sep=':'), snp)
        coef= summary( cox_coef)$coefficients[2]
        sd= summary(cox_coef)$table[2,2]
        n= summary(cox_coef)$n
        pvalue= summary(cox_coef)$table[2,4]
}
        txt = sprintf( "%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\n", snp, ref, eff, 'additive', n, het, hom, coef, sd, pvalue)
cat(txt, file= outfile, append= T)


}
)
}



con= file(infile)
open(con)


repeat {
block.text= readLines(con, n= 250)

if (length(block.text) == 0){ # if there's nothing left, leave the loop
        break
}

funk(block.text)

}

close(con)
}

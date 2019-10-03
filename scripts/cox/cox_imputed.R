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

if (file.size(infile)==0) {
d= data.frame()
write(d, outfile)
quit(status=0)
}


if (grepl('harvest', phenofile)){
pheno= mutate(pheno, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 & 
                INDUKSJON_AMNIOTOMI==0), 
                PARITY0= as.numeric(PARITET_5==0))
} else if (!grepl('harvest', phenofile)){
pheno= mutate(pheno, spont= as.numeric(FSTART=='Spontan' | FSTART== '' & (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'))
names(pheno)[names(pheno) == 'SentrixID'] <- 'SentrixID_1'
}

pheno= select(pheno, SentrixID_1, SVLEN_UL_DG, spont, PARITY0, PC1, PC2, PC3, PC4, PC5, PC6)

ids= readLines(snakemake@input[[3]])

colnames= c(list('chr', 'pos', 'ref', 'eff'), ids)
colnames= unlist(colnames, use.names=FALSE)


funk= function(block.text){
	dataChunk= fread(text=block.text, sep="\t", col.names= colnames)
	if (ncol(dataChunk)==0) return(NULL)
	df= as.matrix(dataChunk[, 5:ncol(dataChunk)], ncol= length(5:ncol(dataChunk)), nrow= nrow(dataChunk))
	df[df== '1|0']= 1
	df[df== '0|1']= 1
	df[df== '0|0']= 0
	df[df== '1|1']= 2

	df= matrix(as.numeric(df), ncol= length(5:ncol(dataChunk)), nrow= nrow(dataChunk))
	eaf= rowSums(df) / (rowSums(!is.na(df)) * 2)
	
	df[which(eaf>0.5), ]= abs(df[which(eaf>0.5), ] -2)
	df[df== 1]= 0
	df[df== 2]= 1
	dataChunk= cbind(dataChunk[,1:4], df)
	dataChunk[which(eaf>0.5) + 4, c(3,4)]= dataChunk[which(eaf>0.5) + 4, c(4,3)]
	names(dataChunk)= colnames
	genvars= paste(dataChunk$chr, dataChunk$pos, dataChunk$ref, dataChunk$eff, sep=':')
	dataChunk= subset(dataChunk, select = -c(chr, pos, ref, eff))

	dataChunk= as.data.frame(t(dataChunk))
	names(dataChunk)[1:length(genvars)]= genvars
	dataChunk= dataChunk[,!apply(dataChunk, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
	if (!is.data.frame(dataChunk)) return(NULL)
        dataChunk$id= gsub('X', '',rownames(dataChunk)) 
        geno= inner_join(pheno, dataChunk, by= c('SentrixID_1' = 'id'))
        cox_coef= lapply(names(geno[,-c(1:dim(pheno)[2])]), function(snp){
        cox_coef= tryCatch(coxph(Surv( geno$SVLEN_UL_DG, geno$spont)~ geno[,snp] + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6, na.action = na.omit), warning = function(cond) {NA}, error = function(cond) {NA})
        if (is.na(cox_coef)) {
        coef= NA
        sd= NA
        n= nrow(df)
        pvalue= NA
	correlation= NA
        corr_pvalue= NA
} else {
        coef = summary( cox_coef)$coefficients[1,1]
        sd= summary(cox_coef)$coefficient[1,3]
        n= summary(cox_coef)$n
        pvalue= summary(cox_coef)$coefficient[1,5]
        zph= cox.zph(cox_coef)
        correlation= zph$table[1, 1]
        corr_pvalue= zph$table[1, 3]
}
	txt = sprintf( "%s\t%e\t%e\t%e\t%e\t%e\t%e\n", snp, n, coef, sd, pvalue, correlation, corr_pvalue)
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

library(dplyr)
library(tidyr)
library(survival)
library(parallel)
library(data.table)

infile= snakemake@input[[1]]
phenofile= snakemake@input[[2]]
ID= 'SentrixID_1' #'MOR_PID' # ID name
outfile= snakemake@output[[1]]
time_t= 'SVLEN_UL_DG' # time variable name
outcome= 'spont' # event variable name
covars= c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PARITY0')

options(stringsAsFactors=FALSE)

pheno= fread(phenofile)

if (grepl('harvest', phenofile)){
pheno= mutate(pheno, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 & 
		INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 & 
		INDUKSJON_AMNIOTOMI==0), 
		PARITY0= as.numeric(PARITET_5==0))
} else if (!grepl('harvest', phenofile)){
pheno= mutate(pheno, spont= as.numeric((FSTART=='Spontan' | FSTART== '') & (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' & 
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' & 
                INDUKSJON_AMNIOTOMI=='Nei'), 
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'))

names(pheno)[names(pheno) == 'SentrixID'] <- 'SentrixID_1'
}

pheno= pheno[!duplicated(pheno$SentrixID_1), ]

# read pheno file; each row is one participant, and column represents one variable
print(infile)

colnames= readLines(gzfile(infile), n= 1)
colnames= unlist(strsplit(colnames, '\t'))


dataChunk= fread(cmd= paste('gzip -cd', infile), sep="\t", col.names= colnames)
        genvars= paste(dataChunk$CHR, dataChunk$segment, sep=':')
	dataChunk= subset( dataChunk, select = -c(CHR,segment))


        dataChunk= as.data.frame(t(dataChunk))
        dataChunk$id= gsub('X','',rownames(dataChunk))
	names(dataChunk)[1:length(genvars)]= genvars
        geno= inner_join(pheno, dataChunk, by= c('SentrixID_1' = 'id'))
	cox_coef= lapply(names(geno[,-c(1:dim(pheno)[2])]), function(snp){cox_coef= survreg(Surv( geno$SVLEN_UL_DG, geno$spont)~ geno[,snp] + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6, na.action = na.omit, dist= 'weibull')
#	coef = summary( cox_coef)$coefficients[1,1]
#	sd= summary(cox_coef)$coefficient[1,3]
	n= summary(cox_coef)$n
#	pvalue= summary(cox_coef)$coefficient[1,5]
#	zph= cox.zph(cox_coef)
#	correlation= zph$table[1, 1]
#	corr_pvalue= zph$table[1, 3]
       coef = summary( cox_coef)$coefficients[2]
       sd= summary(cox_coef)$table[2,2]
       pvalue= summary(cox_coef)$table[2, 4]
	txt = sprintf( "%s\t%e\t%e\t%e\t%e\n", snp, n, coef, sd, pvalue)
cat(txt, file= outfile, append= T)
}
)

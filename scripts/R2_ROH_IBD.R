library(data.table)
library(dplyr)

pca= fread(snakemake@input[[1]])

fam_ibd= fread(snakemake@input[[2]])
names(fam_ibd)= c('FID', 'IID', 'x1','x2', 'x3','x4')

fam_roh= fread(snakemake@input[[3]])
names(fam_roh)= c('FID', 'IID', 'x1','x2', 'x3','x4')

ibd= fread(snakemake@input[[4]])


if (grepl('harvest', snakemake@input[[1]])) {
names(pca)= c('IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
} else {
names(pca)= c('FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
}


infile= list()
R2= list()



flist= snakemake@input[grep('.hom.indiv', snakemake@input)]

for (f in flist) {
	infile= c(f, infile)
	d= fread(f)
	d= full_join(d, ibd, by= c('IID'= 'Child'))

	d= filter(d, Mother %in% pca$IID,
                Father %in% pca$IID,
                Mother %in% fam_ibd$IID,
                IID %in% fam_roh$IID,
                Father %in% fam_ibd$IID,
                Mother %in% fam_ibd$IID,
		IID %in% pca$IID)
	
	d$cM= ifelse(is.na(d$cM), 0, d$cM)
	d$KB= ifelse(is.na(d$KB), 0, d$KB)
	
	
	r= (cor(d$cM, d$KB, use= 'complete'))^2
	R2= c(r, R2)
}

d= do.call(rbind, Map(data.frame, R2=R2, file=infile))

write.table(d, snakemake@output[[1]], sep= '\t', row.names= FALSE, col.names= TRUE, quote= FALSE)

d= d[order(d$R2, decreasing=T), ]

winner= as.character(d[1, 'file'])

prunning= ifelse(grepl('none', winner), 0, ifelse(grepl('soft', winner), 1, ifelse(grepl('moderate', winner), 2, 3)))

winner= strsplit(winner, '_')
winner= lapply(winner, function(x) gsub('.hom.indiv', '', x))
winner= lapply(unlist(winner), as.numeric)
winner= winner[!is.na(winner)]
winner= c(prunning, winner)

write.table(t(as.data.frame(winner)), snakemake@output[[2]], row.names=F, col.names= F, sep= '\t', quote=F)


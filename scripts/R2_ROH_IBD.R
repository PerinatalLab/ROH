library(data.table)
library(dplyr)

pca= fread(snakemake@input[[1]])
pcal= list(pca$IID)

if (grepl('harvest', snakemake@input[[1]])) {
names(pca)= c('IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
} else {
names(pca)= c('FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
}


ibd= unlist(snakemake@input[grepl('ibd', snakemake@input)])


ibd= fread(ibd)
ibd= ibd %>% select(Child, cM)

fam= fread(snakemake@input[[2]])
names(fam)= c('FID', 'IID', 'x1','x2', 'x3','x4')
fam= fam %>% filter(IID %in% pca$IID) %>% select(IID, x1)
fam12= fread(snakemake@input[[3]])
names(fam12)= c('FID', 'IID', 'x1','x2', 'x3','x4')

trio= fread(snakemake@input[[4]])
trio= trio[,2:4]

infile= list()
R2= list()
NSEG= list()
KBAVG= list()

flist= snakemake@input[grep('.hom.indiv', snakemake@input)]

for (f in flist) {
	infile= c(f, infile)
	d= fread(f)
	d= left_join(d, trio, by= c('IID'= 'Child'))
	d= full_join(d, fam, by= 'IID')
	d= left_join(d, ibd, by= c('IID'= 'Child')) %>% inner_join(., pca, by= c('IID'))
	d= d %>% filter(IID %in% pca$IID,
                       Mother %in% pca$IID,
                       Father %in% pca$IID)

	d$cM= ifelse(is.na(d$cM), 0, d$cM)
	d$KB= ifelse(is.na(d$KB), 0, d$KB)
	d$NSEG= ifelse(is.na(d$NSEG), 0, d$NSEG)
	d$KBAVG= ifelse(is.na(d$KBAVG), 0, d$KBAVG)
	if (grepl('harvest', f)) {
	d$BATCH= ifelse(d$IID %in% fam12$IID, 1, 0)
	d$resid= resid(lm(cM~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+ BATCH, d))
	} else {
	d$resid= resid(lm(cM~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, d))
	}
	d= d %>% filter(IID %in% fam$IID)
	r= (cor(d$cM, d$KB, use= 'complete'))^2
	segs= (cor(d$cM, d$NSEG, use= 'complete'))^2
	kbavg= (cor(d$cM, d$KBAVG, use= 'complete'))^2
	R2= c(r, R2)
	NSEG= c(segs, NSEG)
	KBAVG= c(kbavg, KBAVG)
}

d= do.call(rbind, Map(data.frame, R2=R2, file=infile, NSEG= NSEG, KBAVG= KBAVG))

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


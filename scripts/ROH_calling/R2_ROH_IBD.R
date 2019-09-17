library(data.table)
library(dplyr)

SelectRelated= function(kin, sample_list){
 kin= kin %>% filter(KINSHIP>0.0884)
 kin= kin %>% filter(ID1 %in% sample_list, ID2 %in% sample_list)
if (nrow(kin) > 0){
 kin= kin %>% select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2", "ID1", "KINSHIP")
  kin_temp= rbind(kin_temp, kin)
  kin_temp= kin_temp %>% add_count(ID1)
  kin_temp= kin_temp %>% add_count(ID2)
  kin_temp= arrange(kin_temp, n, nn)
  to_keep= list()

  for (i in 1:nrow(kin_temp)) {
    if (kin_temp[i,"ID1"] %in% unlist(kin_temp[0:i,"ID2"])) {
      kin_temp[i,"ID2"]= "X"
    }
    else
      to_keep[[i]] <- kin_temp[["ID1"]][i]
  }
  to_remove= kin_temp %>% filter(!(ID1 %in% unlist(to_keep))) %>% select(ID1)
  to_remove= to_remove[!duplicated(to_remove$ID1),] 

  return(unlist(to_remove[,1]))
}
}


pca= fread(snakemake@input[[1]])

fam_ibd= fread(snakemake@input[[2]])
names(fam_ibd)= c('FID', 'IID', 'x1','x2', 'x3','x4')

fam_roh= fread(snakemake@input[[3]])
names(fam_roh)= c('FID', 'IID', 'x1','x2', 'x3','x4')

ibd= fread(snakemake@input[[4]])

trio= fread(snakemake@input[[6]])
trio= filter(trio, Child %in% fam_ibd$IID, Father %in% fam_ibd$IID, Mother %in% fam_ibd$IID)

ibd= full_join(ibd, trio, by= c('Child', 'Mother', 'Father'))

flag= fread(snakemake@input[[5]])

pca_out= fread(snakemake@input[[8]], h=F)

if (grepl('harvest', snakemake@input[[1]])) {
names(pca)= c('FID', 'IID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
flag= rename(flag, coreLMM = coreOK, phenoOK= phenotypesOK)
} else {
names(pca)= c('FID', 'IID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
}


flag= filter(flag, coreLMM== TRUE, genotypesOK== TRUE, phenoOK== TRUE)

kin= fread(snakemake@input[[7]])

trio= trio %>% filter(Mother %in% flag$IID, Father %in% flag$IID, Child %in% flag$IID)

flag= flag %>% filter(!IID %in% SelectRelated(kin, trio$Mother))
trio= trio %>% filter(Mother %in% flag$IID)
flag= flag %>% filter(!IID %in% SelectRelated(kin, trio$Father))
trio= trio %>% filter(Father %in% flag$IID)
flag= flag %>% filter(!IID %in% SelectRelated(kin, trio$Child))
trio= trio %>% filter(Child %in% flag$IID)

infile= list()
R2= list()



flist= snakemake@input[grep('.hom.indiv', snakemake@input)]

for (f in flist) {
	infile= c(f, infile)
	d= fread(f)
	d= full_join(d, fam_roh, by= c('IID'))
	d= inner_join(d, ibd, by= c('IID'= 'Child'))

	d= filter(d, Mother %in% flag$IID,
                Father %in% flag$IID, 
		IID %in% flag$IID,
		!(IID %in% pca_out$V2),
		!(Father %in% pca_out$V2),
		!(Mother %in% pca_out$V2))
	
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

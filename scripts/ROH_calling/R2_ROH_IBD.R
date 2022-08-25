library(data.table)
library(dplyr)
library(caret)

SelectRelated= function(kin, sample_list){
 kin= kin %>% filter(KINSHIP>0.0884)
 kin= kin %>% filter(ID1 %in% sample_list, ID2 %in% sample_list)
if (nrow(kin) > 0){
 kin= kin %>% dplyr::select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2", "ID1", "KINSHIP")
  kin_temp= rbind(kin_temp, kin)
  kin_temp= kin_temp %>% add_count(ID1, name= 'n_ID1')
  kin_temp= kin_temp %>% add_count(ID2, name= 'n_ID2')
  kin_temp= arrange(kin_temp, n_ID1, n_ID2)
  to_keep= list()

  for (i in 1:nrow(kin_temp)) {
    if (kin_temp[i,"ID1"] %in% unlist(kin_temp[0:i,"ID2"])) {
      kin_temp[i,"ID2"]= "X"
    }
    else
      to_keep[[i]] <- kin_temp[["ID1"]][i]
  }
  to_remove= kin_temp %>% filter(!(ID1 %in% unlist(to_keep))) %>% dplyr::select(ID1)
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
trio= filter(trio, Child %in% fam_roh$IID, Father %in% fam_ibd$IID, Mother %in% fam_ibd$IID)

ibd= full_join(ibd, trio, by= c('Child', 'Mother', 'Father'))

flag= fread(snakemake@input[[5]])

mfr= fread(snakemake@input[[9]])

link= fread(snakemake@input[[10]], h=T)

pca_out= fread(snakemake@input[[8]], h=F)

if (grepl('harvest', snakemake@input[[1]])) {
names(pca)= c('FID', 'IID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
flag= rename(flag, coreLMM = coreOK, phenoOK= phenotypesOK)
link= dplyr::select(link, PREG_ID_1724, SentrixID_1)
mfr= inner_join(mfr, link, by= 'PREG_ID_1724')

mfr= mutate(mfr, FLERFODSEL==0 , DODKAT<6 | DODKAT>10, !is.na(SVLEN_UL_DG), SVLEN_UL_DG<308, SVLEN_UL_DG>154, is.na(IVF), FOSTERV_POLYHYDRAMNION==0, C00_MALF_ALL==0, FOSTERV_OLIGOHYDRAMNION== 0, VEKT>1500, !is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK==0, HYPERTENSJON_ALENE==0, !is.na(PREEKL))

flag= filter(flag, IID %in% mfr$SentrixID_1)


} else {
names(pca)= c('FID', 'IID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
print(names(link))
link= dplyr::select(link, PREG_ID_315, SentrixID)

mfr= inner_join(mfr, link, by= 'PREG_ID_315')

mfr= mutate(mfr, FLERFODSEL== 'Enkeltfødsel' , grepl('Levendefødt', DODKAT), !is.na(SVLEN_UL_DG),SVLEN_UL_DG<308, SVLEN_UL_DG>154, is.na(IVF) | IVF== '', FOSTERV_POLYHYDRAMNION=='Nei', C00_MALF_ALL=='Nei', FOSTERV_OLIGOHYDRAMNION== 'Nei', VEKT>1500, !is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK==0, HYPERTENSJON_ALENE==0, !is.na(PREEKL))

flag= filter(flag, IID %in% mfr$SentrixID)

}


flag= filter(flag, genotypesOK== TRUE, phenoOK== TRUE)

kin= fread(snakemake@input[[7]])
ibd= ibd %>% filter(Mother %in% flag$IID, Father %in% flag$IID, Child %in% flag$IID)

infile= list()
R2= list()
aic_list= list()
out_list= snakemake@input[grep('list_', snakemake@input)]
flist= snakemake@input[grep('.hom.indiv', snakemake@input)]

for (f in flist) {
	for (outUPD in out_list){
	if (sub('.txt', '.hom.indiv', sub('list_', '', outUPD)) == f){
	x= readLines(outUPD)
	}
	}
	infile= c(f, infile)
	d= fread(f)
	d= filter(d, !(IID %in% x))
	d= full_join(d, ibd, by= c('IID'= 'Child'))

	d= filter(d, Mother %in% flag$IID,
                Father %in% flag$IID, 
		IID %in% flag$IID,
		!(IID %in% pca_out$V2),
		!(Father %in% pca_out$V2),
		!(Mother %in% pca_out$V2))
	d= filter(d, !Mother %in% SelectRelated(kin, d$Mother))
	d= filter(d, !Father %in% SelectRelated(kin, d$Father))
	d= filter(d, !IID %in% SelectRelated(kin, d$IID))
	
	d$cM= ifelse(is.na(d$cM), 0, d$cM)
	d$KB= ifelse(is.na(d$KB), 0, d$KB)
	d= inner_join(d, pca, by= 'IID')
#	set.seed(0)
#	train_control <- trainControl(method="repeatedcv", number=10, repeats= 100)
#	d$rescM= lm(cM~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, d)$resid
#	model <- train(rescM~ KB, data=d, trControl=train_control, method="lm", na.action= na.omit)
#	r= summary(lm(cM~ KB, data=d, na.action= na.omit))$r.squared
#	r= mean(model$resample$Rsquared, na.rm= T)
	r= with(d, cor(d$cM, d$KB, use= 'complete', method= 'spearman'))
	R2[[f]]= r
}
d= data.frame(do.call('rbind', R2))
d$file= rownames(d)
names(d)= c('R2', 'file')

write.table(d, snakemake@output[[1]], sep= '\t', row.names= F, col.names= T, quote= FALSE)

d= d[order(d$R2, decreasing=T), ]

winner= as.character(d[1, 'file'])

prunning= ifelse(grepl('none', winner), 0, ifelse(grepl('soft', winner), 1, ifelse(grepl('moderate', winner), 2, 3)))

winner= strsplit(winner, '_')
winner= lapply(winner, function(x) gsub('.hom.indiv', '', x))
winner= lapply(unlist(winner), as.numeric)
winner= winner[!is.na(winner)]
winner= c(prunning, winner)

write.table(t(as.data.frame(winner)), snakemake@output[[2]], row.names=F, col.names= F, sep= '\t', quote=F)


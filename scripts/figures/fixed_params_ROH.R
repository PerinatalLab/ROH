library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(caret)

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')

SelectRelated= function(kin, sample_list){
 kin= kin %>% filter(KINSHIP>0.0884)
 kin= kin %>% filter(ID1 %in% sample_list, ID2 %in% sample_list)
if (nrow(kin) > 0){
 kin= kin %>% select(ID1, ID2, KINSHIP)
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
  to_remove= kin_temp %>% filter(!(ID1 %in% unlist(to_keep))) %>% select(ID1)
  to_remove= to_remove[!duplicated(to_remove$ID1),]

  return(unlist(to_remove[,1]))
}
}

input= unlist(snakemake@input)
cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

df_list= list()
r_list= list()

for (coh in cohorts){


input_coh= input[grep(coh, input)]

out_list= input_coh[grep('list_', input_coh)]

pca= fread(unlist(input_coh[grep('pca.txt', input_coh)]))

fam_ibd= fread(unlist(input_coh[grepl('ibd/to_phase.fam', input_coh)]))
names(fam_ibd)= c('FID', 'IID', 'x1','x2', 'x3','x4')

fam_roh= fread(unlist(input_coh[grepl('fetal.fam', input_coh)]))
names(fam_roh)= c('FID', 'IID', 'x1','x2', 'x3','x4')

ibd= fread(unlist(input_coh[grepl('parental_ibd.txt', input_coh)]))

trio= fread(unlist(input_coh[grepl('trios.txt', input_coh)]))
trio= filter(trio, Child %in% fam_roh$IID, Father %in% fam_ibd$IID, Mother %in% fam_ibd$IID)

ibd= full_join(ibd, trio, by= c('Child', 'Mother', 'Father'))

flag= fread(unlist(input_coh[grepl('flag_list.txt', input_coh)]))

pca_out= fread(unlist(input_coh[grepl('exclude', input_coh)]), h=F)
mfr= fread(unlist(input_coh[grepl('mfr', input_coh)]))
link= fread(unlist(input_coh[grepl('linkage', input_coh)]))

if (coh== 'harvestm12' | coh == 'harvestm24') {
flag= rename(flag, phenoOK= phenotypesOK)
link= select(link, PREG_ID_1724, SentrixID_1)
mfr= inner_join(mfr, link, by= 'PREG_ID_1724')

mfr= mutate(mfr, FLERFODSEL==0 , DODKAT<6 | DODKAT>10, !is.na(SVLEN_UL_DG), SVLEN_UL_DG<308, SVLEN_UL_DG>154, is.na(IVF), FOSTERV_POLYHYDRAMNION==0, C00_MALF_ALL==0, FOSTERV_OLIGOHYDRAMNION== 0, VEKT>1500)

flag= filter(flag, IID %in% mfr$SentrixID_1)

} else {

link= select(link, PREG_ID_315, SentrixID)

mfr= inner_join(mfr, link, by= 'PREG_ID_315')
mfr= mutate(mfr, FLERFODSEL== 'Enkeltfødsel' , grepl('Levendefødt', DODKAT), !is.na(SVLEN_UL_DG),SVLEN_UL_DG<308, SVLEN_UL_DG>154, is.na(IVF) | IVF== '', FOSTERV_POLYHYDRAMNION=='Nei', C00_MALF_ALL=='Nei', FOSTERV_OLIGOHYDRAMNION== 'Nei', VEKT>1500)

flag= filter(flag, IID %in% mfr$SentrixID)

}

names(pca)= c('FID', 'IID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

pca= select(pca, IID, PC1, PC2, PC3, PC4, PC5, PC6)

flag= filter(flag, genotypesOK== TRUE, phenoOK== TRUE)

kin= fread(unlist(input_coh[grepl('.kin0', input_coh)]))

trio= trio %>% filter(Mother %in% flag$IID, Father %in% flag$IID, Child %in% flag$IID)

out_list= input_coh[grep('list_', input_coh)]


x= readLines(unlist(input_coh[grepl('list_', input_coh)]))
d= fread(unlist(input_coh[grepl('hom.indiv', input_coh)]))
d= left_join(d, pca, by= 'IID')
d= full_join(d, fam_roh, by= c('IID'))
d= inner_join(d, ibd, by= c('IID'= 'Child'))
d= filter(d, !(IID %in% x))
d= filter(d, Mother %in% flag$IID,
        Father %in% flag$IID,
        IID %in% flag$IID,
	!(IID %in% x),
        !(IID %in% pca_out$V2),
        !(Father %in% pca_out$V2),
        !(Mother %in% pca_out$V2))
	d= filter(d, !Mother %in% SelectRelated(kin, d$Mother))
        d= filter(d, !Father %in% SelectRelated(kin, d$Father))
        d= filter(d, !IID %in% SelectRelated(kin, d$IID))

d$cM= ifelse(is.na(d$cM), 0, d$cM)
d$KB= ifelse(is.na(d$KB), 0, d$KB)

d= select(d, cM, KB, PC1, PC2, PC3, PC4, PC5, PC6)
d$cohort= coh

train_control <- trainControl(method="cv", number=10)
d$rescM= lm(cM~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, d)$resid
model <- train(rescM~ KB, data=d, trControl=train_control, method="lm", na.action= na.omit)

r_list[[coh]]= mean(model$resample$Rsquared)

df_list= c(df_list, list(d))
}




d= do.call(bind_rows, df_list)

d$cohort2= ifelse((d$cohort== 'harvestm12' | d$cohort== 'harvestm24' | d$cohort== 'rotterdam1'), 0, 1)

d$cohort= factor(d$cohort, levels= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay'))

#d$cM= ifelse(d$cM==0, 10**-6, d$cM)



lab= as.data.frame(do.call('rbind', r_list))
lab$cohort= rownames(lab)
names(lab)= c('R', 'cohort')


print(group_by(d, cohort) %>% summarize(n= sum(!is.na(cM) & !is.na(KB))))
print(lab)

x1= ggplot(data= d, aes(x= cM, y= KB/10**6, colour= cohort)) +
geom_point(size= 1.5) +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name = "Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= c(colors_3, 'black', 'grey', 'cyan')) +
#scale_shape_manual(name="Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= rep(15:16, 3)) +
geom_smooth(method = "lm", se=FALSE, colour= "black", formula = y ~ x, size= 0.6, linetype = 'dashed') +
facet_wrap(vars(cohort), ncol= 3) +
geom_text(data=lab, aes(x= Inf, y= -Inf, label= paste0('R**2:  ', sprintf('%0.2f', round(R,2)))), hjust= 1, vjust= -1, parse= T, colour= 'black', show.legend = FALSE) +
theme(strip.text = element_blank(),
          strip.background = element_blank()) +
xlab('Parental genetic relatedness, total cM') +
    ylab('Offspring ROH length, cM')

save_plot(snakemake@output[[1]], plot= x1, base_width=297, base_height=210, units="mm", device= cairo_ps)


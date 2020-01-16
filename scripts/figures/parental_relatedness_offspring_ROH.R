library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')

SelectRelated= function(kin, sample_list){
 kin= kin %>% filter(KINSHIP>0.0884)
 kin= kin %>% filter(ID1 %in% sample_list, ID2 %in% sample_list)
if (nrow(kin) > 0){
 kin= kin %>% select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2", "ID1", "KINSHIP")
  kin_temp= rbind(kin_temp, kin)
  kin_temp= kin_temp %>% add_count(ID1, name= 'n')
  kin_temp= kin_temp %>% add_count(ID2, name= 'nn')
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

input= unlist(snakemake@input)
cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

df_list= list()

for (coh in cohorts){


input_coh= input[grep(coh, input)]

pca= fread(input_coh[grep('pca.txt', input_coh)])

fam_ibd= fread(unlist(input_coh[grepl('ibd/to_phase.fam', input_coh)]))
names(fam_ibd)= c('FID', 'IID', 'x1','x2', 'x3','x4')

fam_roh= fread(unlist(input_coh[grepl('fetal.fam', input_coh)]))
names(fam_roh)= c('FID', 'IID', 'x1','x2', 'x3','x4')

ibd= fread(unlist(input_coh[grepl('parental_ibd.txt', input_coh)]))

trio= fread(unlist(input_coh[grepl('trios.txt', input_coh)]))
trio= filter(trio, Child %in% fam_ibd$IID, Father %in% fam_ibd$IID, Mother %in% fam_ibd$IID)

ibd= full_join(ibd, trio, by= c('Child', 'Mother', 'Father'))

flag= fread(unlist(input_coh[grepl('flag_list.txt', input_coh)]))

pca_out= fread(unlist(input_coh[grepl('all_pca_outliers_hapmap.txt', input_coh)]), h=F)

if (coh== 'harvestm12' | coh == 'harvestm24') {
flag= rename(flag, coreLMM = coreOK, phenoOK= phenotypesOK)
}

names(pca)= c('FID', 'IID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

flag= filter(flag, coreLMM== TRUE, genotypesOK== TRUE, phenoOK== TRUE)

kin= fread(unlist(input_coh[grepl('.kin0', input_coh)]))

trio= trio %>% filter(Mother %in% flag$IID, Father %in% flag$IID, Child %in% flag$IID)

flag= flag %>% filter(!IID %in% SelectRelated(kin, trio$Mother))
trio= trio %>% filter(Mother %in% flag$IID)
flag= flag %>% filter(!IID %in% SelectRelated(kin, trio$Father))
trio= trio %>% filter(Father %in% flag$IID)
flag= flag %>% filter(!IID %in% SelectRelated(kin, trio$Child))
trio= trio %>% filter(Child %in% flag$IID)


arg= fread(unlist(input_coh[grepl('arg_R2', input_coh)]))
arg= arg[arg$R2== max(arg$R2, na.rm=T), ]
d= fread(arg$file)

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

d= select(d, cM, KB)
d$cohort= coh

df_list= c(df_list, list(d))
}




d= do.call(bind_rows, df_list)

d$cohort2= ifelse((d$cohort== 'harvestm12' | d$cohort== 'harvestm24' | d$cohort== 'rotterdam1'), 0, 1)

d$cohort= factor(d$cohort, levels= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay'))

r_list= c()
coh_list= c()
for (coh in cohorts){
r_list= c(r_list, with(d[d$cohort== coh, ], cor(cM, KB, use= 'complete')**2))
coh_list= c(coh_list, coh)
}



lab= data.frame(R= r_list, cohort= coh_list)

print(group_by(d, cohort) %>% summarize(n= sum(!is.na(cM) & !is.na(KB))))

print(cor(d$cM, d$KB, use= 'complete'))

x1= ggplot(d, aes(x= cM, y= KB)) +
geom_point(aes( colour= cohort, shape= cohort), size = 1.5) + 
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name = "Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= rep(colors_3, 2)) +
scale_shape_manual(name="Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= rep(15:16, 3)) +
geom_smooth(method = "lm", se=FALSE, colour= "black", formula = y ~ x, size= 0.6, linetype = 'dashed') +
facet_wrap(vars(cohort), ncol= 3) +
geom_text(data=lab, aes(x= Inf, y= -Inf, label= paste0('R**2:  ', sprintf('%0.2f', round(R,2)))), hjust= 1, vjust= -1, parse= T) +
theme(strip.text = element_blank(),
          strip.background = element_blank()) +
xlab('Parental genetic relatedness, total cM') +
    ylab('Offspring ROH length, cM')

save_plot(snakemake@output[[1]], plot= x1, base_width=297, base_height=210, units="mm")






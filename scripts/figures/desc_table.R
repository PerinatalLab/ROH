library(dplyr)
library(data.table)
library(metafor)
library(survival)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

mom_list= list()

input= unlist(snakemake@input)

for (coh in cohorts){
input_coh= input[grep(coh, input)]

mom= fread(input_coh[grep('mfr_maternal', input_coh)])
q1= fread(input_coh[grep('q1', input_coh)])

if (coh== 'harvestm12' | coh== 'harvestm24'){
mom= rename(mom, 'PREG_ID_315'= 'PREG_ID_1724')
q1= rename(q1, 'PREG_ID_315'= 'PREG_ID_1724')

mom= mutate(mom, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 &
                INDUKSJON_AMNIOTOMI==0),
                PARITY0= as.numeric(PARITET_5==0))
}

if (coh!= 'harvestm12' & coh!= 'harvestm24'){
mom= mutate(mom, spont= as.numeric((FSTART=='Spontan' | FSTART== '') & ((KSNITT=='') | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'))

q1$AA1124= as.integer(factor(q1$AA1124, c("9-årig grunnskole", "1-2-årig videregående", "Videregående yrkesfaglig", "3-årig videregående allmennfaglig, gymnas", "Distriktshøyskole, universitet inntil 4 år (cand. mag., sykepleier, lærer, ingeniør)", "Universitet, høyskole, mer enn 4 år (hovedfag, embetseksamen)")))
q1$AA1125= as.integer(factor(q1$AA1125, c("9-årig grunnskole", "1-2-årig videregående", "Videregående yrkesfaglig", "3-årig videregående allmennfaglig, gymnas", "Distriktshøyskole, universitet inntil 4 år (cand. mag., sykepleier, lærer, ingeniør)", "Universitet, høyskole, mer enn 4 år (hovedfag, embetseksamen)")))
q1$AA1315= as.integer(factor(q1$AA1315, c("Ingen Inntekt", "Under 150.000 kr.", "150-199.999 kr.", "200-299.999 kr.", "300-399.999 kr.", "400-499.999 kr.", "over 500.000 kr.")))
q1$AA1316= as.integer(factor(q1$AA1316, c("Ingen Inntekt", "Under 150.000 kr.", "150-199.999 kr.", "200-299.999 kr.", "300-399.999 kr.", "400-499.999 kr.", "over 500.000 kr.")))
mom$MORS_ALDER= mom$FAAR - mom$MOR_FAAR
}

q1$edu[q1$AA1124 <3 | q1$AA1125<3]= 1
q1$edu[q1$AA1124 ==3 | q1$AA1125==3]= 2
q1$edu[q1$AA1124 ==4 | q1$AA1125==4]= 2
q1$edu[q1$AA1124 ==5 | q1$AA1125==5]= 3
q1$edu[q1$AA1124 ==6 | q1$AA1125==6]= 4
q1$edu[(is.na(q1$edu) & q1$AA1128==6) | (is.na(q1$edu) & q1$AA1129==1)]= 5

q1$AA1315[q1$AA1315==1 | q1$AA1315==2| q1$AA1315==3| q1$AA1315==4]= 0
q1$AA1315[q1$AA1315==5 | q1$AA1315==6 | q1$AA1315==7]= 1
q1$AA1316[q1$AA1316==1 | q1$AA1316==2| q1$AA1316==3| q1$AA1316==4]= 0
q1$AA1316[q1$AA1316==5 | q1$AA1316==6 | q1$AA1316==7]= 1

q1$income[q1$AA1315==0 & q1$AA1316==0]= 0
q1$income[q1$AA1315==0 & is.na(q1$AA1316)]= 0
q1$income[is.na(q1$AA1315) & q1$AA1316==0]= 0
q1$income[q1$AA1315==1 | q1$AA1316==1]= 1
q1$income[q1$AA1315==1 & is.na(q1$AA1316)]= 1
q1$income[is.na(q1$AA1315) & q1$AA1316==1]= 1
q1$income[q1$AA1315==1 & q1$AA1316==1]= 2

q1$edu= ifelse(q1$edu>=4, 4, q1$edu)
q1$edu= ifelse(q1$edu<2, 2, q1$edu)
q1$edu= q1$edu -1 

q1= select(q1, PREG_ID_315, income, edu)

mom= left_join(mom, q1, by= 'PREG_ID_315')
mom= select(mom, spont, SVLEN_UL_DG, PARITY0, VEKT, MORS_ALDER, income, edu)
mom$cohort= coh

mom_list= c(mom_list, list(mom))

}

mom= do.call('rbind', mom_list)


newlist= lapply(cohorts, function(x) {
df= mom[mom$cohort== x, ]
x_list= data.frame(
cohort= x,
total_n= length(df$spont),
spont_n= sum(df$spont, na.rm= T),
spont_frac= round(mean(df$spont, na.rm= T), 2),
ga_mean= round(mean(df$SVLEN_UL_DG, na.rm= T)),
ga_sd= round(sd(df$SVLEN_UL_DG, na.rm= T), 1),
parity_n= sum(df$PARITY0, na.rm=T),
parity_frac= round(mean(df$PARITY0, na.rm= T), 2),
bw_mean= round(mean(df$VEKT, na.rm= T)),
bw_sd= round(sd(df$VEKT, na.rm= T)),
age_mean= round(mean(df$MORS_ALDER, na.rm= T)),
age_sd= round(sd(df$MORS_ALDER, na.rm= T), 1),
inc1_n= sum(df$income== 0, na.rm= T),
inc1_frac= round(mean(df$income== 0, na.rm= T), 2),
inc2_n= sum(df$income== 1, na.rm= T),
inc2_frac= round(mean(df$income== 1, na.rm= T), 2),
inc3_n= sum(df$income== 2, na.rm= T),
inc3_frac= round(mean(df$income== 2, na.rm= T) , 2),
edu1_n= sum(df$edu== 1, na.rm= T),
edu1_frac= round(mean(df$edu== 1, na.rm= T), 2),
edu2_n= sum(df$edu== 2, na.rm= T),
edu2_frac= round(mean(df$edu== 2, na.rm= T), 2),
edu3_n= sum(df$edu== 3, na.rm= T),
edu3_frac= round(mean(df$edu== 3, na.rm= T),2))
return(x_list)
})

x= do.call('rbind', newlist)

spont_pvalue= chisq.test(table(mom$spont, mom$cohort))$p.value
parity_pvalue= chisq.test(table(mom$PARITY0, mom$cohort))$p.value
ga_pvalue= summary(aov(SVLEN_UL_DG~ cohort, mom))[[1]][["Pr(>F)"]][1]
bw_pvalue= summary(aov(VEKT~ cohort, mom))[[1]][["Pr(>F)"]][1]
age_pvalue= summary(aov(MORS_ALDER~ cohort, mom))[[1]][["Pr(>F)"]][1]
inc_pvalue= chisq.test(table(mom$income, mom$cohort))$p.value
edu_pvalue= chisq.test(table(mom$edu, mom$cohort))$p.value

z= data.frame(var= c('spont', 'ga', 'parity', 'bw', 'age', 'income', 'edu'), pvalue= c(spont_pvalue, ga_pvalue, parity_pvalue, bw_pvalue, age_pvalue, inc_pvalue, edu_pvalue))


write.table(x, snakemake@output[[1]], sep= '\t', row.names=F, col.names=T, quote=F)
write.table(z, snakemake@output[[2]], sep= '\t', row.names=F, col.names=T, quote=F)


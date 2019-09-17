library(dplyr)
library(data.table)
library(metafor)
library(survival)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

mom_list= list()
dad_list= list()
fet_list= list()

input= unlist(snakemake@input)

for (coh in cohorts){
input_coh= unlist(input[grep(coh, input)])

q1= fread(unlist(input_coh[grep('q1', input_coh)]))

mom= fread(unlist(input_coh[grep('mfr_maternal', input_coh)]))
dad= fread(unlist(input_coh[grep('mfr_paternal', input_coh)]))
fet= fread(unlist(input_coh[grep('mfr_fetal', input_coh)]))


if (coh== 'harvestm12' | coh== 'harvestm24'){
mom= rename(mom, 'PREG_ID_315'= 'PREG_ID_1724')
dad= rename(dad, 'PREG_ID_315'= 'PREG_ID_1724')
fet= rename(fet, 'PREG_ID_315'= 'PREG_ID_1724')
q1= rename(q1, 'PREG_ID_315'= 'PREG_ID_1724')


mom= mutate(mom, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 &
                INDUKSJON_AMNIOTOMI==0),
                PARITY0= as.numeric(PARITET_5==0))
dad= mutate(dad, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 &
                INDUKSJON_AMNIOTOMI==0),
                PARITY0= as.numeric(PARITET_5==0))
fet= mutate(fet, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 &
                INDUKSJON_AMNIOTOMI==0),
                PARITY0= as.numeric(PARITET_5==0))
}

if (coh!= 'harvestm12' & coh!= 'harvestm24'){
mom= mutate(mom, spont= as.numeric(FSTART=='Spontan' & ((KSNITT=='') | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'))
dad= mutate(dad, spont= as.numeric(FSTART=='Spontan' & ((KSNITT=='') | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'))
fet= mutate(fet, spont= as.numeric(FSTART=='Spontan' & ((KSNITT=='') | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'))

q1$AA1124= as.integer(factor(q1$AA1124, c("9-årig grunnskole", "1-2-årig videregående", "Videregående yrkesfaglig", "3-årig videregående allmennfaglig, gymnas", "Distriktshøyskole, universitet inntil 4 år (cand. mag., sykepleier, lærer, ingeniør)", "Universitet, høyskole, mer enn 4 år (hovedfag, embetseksamen)")))
q1$AA1125= as.integer(factor(q1$AA1125, c("9-årig grunnskole", "1-2-årig videregående", "Videregående yrkesfaglig", "3-årig videregående allmennfaglig, gymnas", "Distriktshøyskole, universitet inntil 4 år (cand. mag., sykepleier, lærer, ingeniør)", "Universitet, høyskole, mer enn 4 år (hovedfag, embetseksamen)")))
q1$AA1315= as.integer(factor(q1$AA1315, c("Ingen Inntekt", "Under 150.000 kr.", "150-199.999 kr.", "200-299.999 kr.", "300-399.999 kr.", "400-499.999 kr.", "over 500.000 kr.")))
q1$AA1316= as.integer(factor(q1$AA1316, c("Ingen Inntekt", "Under 150.000 kr.", "150-199.999 kr.", "200-299.999 kr.", "300-399.999 kr.", "400-499.999 kr.", "over 500.000 kr.")))
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

q1= select(q1, PREG_ID_315, income, edu)

mom= left_join(mom, q1, by= 'PREG_ID_315')
dad= left_join(dad, q1, by= 'PREG_ID_315')
fet= left_join(fet, q1, by= 'PREG_ID_315')

d_mom= fread(unlist(input_coh[grep('maternal_excess', input_coh)]))
d_dad= fread(unlist(input_coh[grep('paternal_excess', input_coh)]))
d_fet= fread(unlist(input_coh[grep('fetal_excess', input_coh)]))

mom= inner_join(mom, d_mom, by= 'IID')
dad= inner_join(dad, d_dad, by= 'IID')
fet= inner_join(fet, d_fet, by= 'IID')

mom= select(mom, none_F, PC1, PC2, PC3, PC4, PC5, PC6, spont, SVLEN_UL_DG, PARITY0, income, edu)
dad= select(dad, none_F, PC1, PC2, PC3, PC4, PC5, PC6, spont, SVLEN_UL_DG, PARITY0, income, edu)
fet= select(fet, none_F, PC1, PC2, PC3, PC4, PC5, PC6, spont, SVLEN_UL_DG, PARITY0, income, edu)


mom$cohort= coh
dad$cohort= coh
fet$cohort= coh

mom_list= c(mom_list, list(mom))
dad_list= c(dad_list, list(dad))
fet_list= c(fet_list, list(fet))

}

mom= do.call('rbind', mom_list)
dad= do.call('rbind', dad_list)
fet= do.call('rbind', fet_list)

df_l= lapply(setNames(cohorts, cohorts), function(coh){
x= coxph(Surv(SVLEN_UL_DG, spont)~ none_F, mom %>% filter(cohort== coh), na.action= na.omit)
return(c(summary(x)$coefficients[1,], n= x$n, event= x$nevent))})
moms_res= as.data.frame(do.call('rbind', df_l))
moms_res$cohort= rownames(moms_res)
names(moms_res)= c('beta', 'HR', 'se', 'z', 'pvalue', 'n', 'event', 'cohort')
moms_res$lowerci = (-1.96* moms_res$se) + moms_res$beta
moms_res$upperci = (1.96* moms_res$se) + moms_res$beta

random= rma(moms_res$beta, sei= moms_res$se)
meta.res= data.frame(cohort= 'Summary effect', beta= random$beta[1], se= random$se[1], pvalue= random$pval[1], lowerci= random$ci.lb, upperci= random$ci.ub, n= sum(moms_res$n), event= sum(moms_res$event))

moms_res= bind_rows(moms_res, meta.res)
moms_res$HR= exp(moms_res$beta)

moms_res$cohort= factor(moms_res$cohort, levels= c('Summary effect', 'normentmay', 'normentfeb', 'rotterdam2', 'rotterdam1', 'harvestm24', 'harvestm12'))

moms_res$sample= 'Mothers'

moms_res$lowerci = exp((-1.96* moms_res$se)+ moms_res$beta)
moms_res$upperci = exp((1.96* moms_res$se)+ moms_res$beta)

df_l= lapply(setNames(cohorts, cohorts), function(coh){
x= coxph(Surv(SVLEN_UL_DG, spont)~  none_F, dad %>% filter(cohort== coh), na.action= na.omit)
return(c(summary(x)$coefficients[1,], n= x$n, event= x$nevent))})
dads_res= as.data.frame(do.call('rbind', df_l))
dads_res$cohort= rownames(dads_res)
names(dads_res)= c('beta', 'HR', 'se', 'z', 'pvalue', 'n', 'event', 'cohort')
dads_res$lowerci = (-1.96* dads_res$se)+ dads_res$beta
dads_res$upperci = (1.96* dads_res$se)+ dads_res$beta

random= rma(dads_res$beta, sei= dads_res$se)
meta.res= data.frame(cohort= 'Summary effect', beta= random$beta[1], se= random$se[1], pvalue= random$pval[1], lowerci= random$ci.lb, upperci= random$ci.ub, n= sum(dads_res$n), event= sum(dads_res$event))

dads_res= bind_rows(dads_res, meta.res)
dads_res$HR= exp(dads_res$beta)
dads_res$cohort= factor(dads_res$cohort, levels= c('Summary effect', 'normentmay', 'normentfeb', 'rotterdam2', 'rotterdam1', 'harvestm24', 'harvestm12'))


dads_res$lowerci = exp((-1.96* dads_res$se)+ dads_res$beta)
dads_res$upperci = exp((1.96* dads_res$se)+ dads_res$beta)

dads_res$sample= 'Fathers'

df_l= lapply(setNames(cohorts, cohorts), function(coh){
x= coxph(Surv(SVLEN_UL_DG, spont)~ none_F, fet %>% filter(cohort== coh), na.action= na.omit)
return(c(summary(x)$coefficients[1,], n= x$n, event= x$nevent))})
fets_res= as.data.frame(do.call('rbind', df_l))
fets_res$cohort= rownames(fets_res)
names(fets_res)= c('beta', 'HR', 'se', 'z', 'pvalue', 'n', 'event', 'cohort')
fets_res$lowerci = (-1.96* fets_res$se)+ fets_res$beta
fets_res$upperci = (1.96* fets_res$se)+ fets_res$beta

random= rma(fets_res$beta, sei= fets_res$se)
meta.res= data.frame(cohort= 'Summary effect', beta= random$beta[1], se= random$se[1], pvalue= random$pval[1], lowerci= random$ci.lb, upperci= random$ci.ub, n= sum(fets_res$n), event= sum(fets_res$event))

fets_res= bind_rows(fets_res, meta.res)
fets_res$HR= exp(fets_res$beta)
fets_res$cohort= factor(fets_res$cohort, levels= c('Summary effect', 'normentmay', 'normentfeb', 'rotterdam2', 'rotterdam1', 'harvestm24', 'harvestm12'))

fets_res$lowerci = exp((-1.96* fets_res$se)+ fets_res$beta)
fets_res$upperci = exp((1.96* fets_res$se)+ fets_res$beta)

fets_res$sample= 'Fetal'

all_res= do.call("rbind", list(moms_res, dads_res, fets_res))

all_res$sample= factor(all_res$sample, levels= c('Mothers','Fathers', 'Fetal'))

write.table(all_res, snakemake@output[[1]], sep= '\t', header= True, quote=F)

mom= do.call('rbind', mom_list)
dad= do.call('rbind', dad_list)
fet= do.call('rbind', fet_list)

df_l= lapply(setNames(cohorts, cohorts), function(coh){
x= coxph(Surv(SVLEN_UL_DG, spont)~ none_F + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY0 + edu + income, mom %>% filter(cohort== coh), na.action= na.omit)
return(c(summary(x)$coefficients[1,], n= x$n, event= x$nevent))})
moms_res= as.data.frame(do.call('rbind', df_l))
moms_res$cohort= rownames(moms_res)
names(moms_res)= c('beta', 'HR', 'se', 'z', 'pvalue', 'n', 'event', 'cohort')
moms_res$lowerci = (-1.96* moms_res$se) + moms_res$beta
moms_res$upperci = (1.96* moms_res$se) + moms_res$beta

random= rma(moms_res$beta, sei= moms_res$se)
meta.res= data.frame(cohort= 'Summary effect', beta= random$beta[1], se= random$se[1], pvalue= random$pval[1], lowerci= random$ci.lb, upperci= random$ci.ub, n= sum(moms_res$n), event= sum(moms_res$event))

moms_res= bind_rows(moms_res, meta.res)
moms_res$HR= exp(moms_res$beta)

moms_res$cohort= factor(moms_res$cohort, levels= c('Summary effect', 'normentmay', 'normentfeb', 'rotterdam2', 'rotterdam1', 'harvestm24', 'harvestm12'))

moms_res$sample= 'Mothers'

moms_res$lowerci = exp((-1.96* moms_res$se)+ moms_res$beta)
moms_res$upperci = exp((1.96* moms_res$se)+ moms_res$beta)

df_l= lapply(setNames(cohorts, cohorts), function(coh){
x= coxph(Surv(SVLEN_UL_DG, spont)~  none_F + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY0 + edu + income, dad %>% filter(cohort== coh), na.action= na.omit)
return(c(summary(x)$coefficients[1,], n= x$n, event= x$nevent))})
dads_res= as.data.frame(do.call('rbind', df_l))
dads_res$cohort= rownames(dads_res)
names(dads_res)= c('beta', 'HR', 'se', 'z', 'pvalue', 'n', 'event', 'cohort')
dads_res$lowerci = (-1.96* dads_res$se)+ dads_res$beta
dads_res$upperci = (1.96* dads_res$se)+ dads_res$beta

random= rma(dads_res$beta, sei= dads_res$se)
meta.res= data.frame(cohort= 'Summary effect', beta= random$beta[1], se= random$se[1], pvalue= random$pval[1], lowerci= random$ci.lb, upperci= random$ci.ub, n= sum(dads_res$n), event= sum(dads_res$event))

dads_res= bind_rows(dads_res, meta.res)
dads_res$HR= exp(dads_res$beta)
dads_res$cohort= factor(dads_res$cohort, levels= c('Summary effect', 'normentmay', 'normentfeb', 'rotterdam2', 'rotterdam1', 'harvestm24', 'harvestm12'))

dads_res$lowerci = exp((-1.96* dads_res$se)+ dads_res$beta)
dads_res$upperci = exp((1.96* dads_res$se)+ dads_res$beta)

dads_res$sample= 'Fathers'

df_l= lapply(setNames(cohorts, cohorts), function(coh){
x= coxph(Surv(SVLEN_UL_DG, spont)~ none_F + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY0 + edu + income, fet %>% filter(cohort== coh), na.action= na.omit)
return(c(summary(x)$coefficients[1,], n= x$n, event= x$nevent))})
fets_res= as.data.frame(do.call('rbind', df_l))
fets_res$cohort= rownames(fets_res)
names(fets_res)= c('beta', 'HR', 'se', 'z', 'pvalue', 'n', 'event', 'cohort')
fets_res$lowerci = (-1.96* fets_res$se)+ fets_res$beta
fets_res$upperci = (1.96* fets_res$se)+ fets_res$beta

random= rma(fets_res$beta, sei= fets_res$se)
meta.res= data.frame(cohort= 'Summary effect', beta= random$beta[1], se= random$se[1], pvalue= random$pval[1], lowerci= random$ci.lb, upperci= random$ci.ub, n= sum(fets_res$n), event= sum(fets_res$event))

fets_res= bind_rows(fets_res, meta.res)
fets_res$HR= exp(fets_res$beta)
fets_res$cohort= factor(fets_res$cohort, levels= c('Summary effect', 'normentmay', 'normentfeb', 'rotterdam2', 'rotterdam1', 'harvestm24', 'harvestm12'))

fets_res$lowerci = exp((-1.96* fets_res$se)+ fets_res$beta)
fets_res$upperci = exp((1.96* fets_res$se)+ fets_res$beta)

fets_res$sample= 'Fetal'

all_res= do.call("rbind", list(moms_res, dads_res, fets_res))

all_res$sample= factor(all_res$sample, levels= c('Mothers','Fathers', 'Fetal'))

write.table(all_res, snakemake@output[[2]], sep= '\t', header= True, quote=F)

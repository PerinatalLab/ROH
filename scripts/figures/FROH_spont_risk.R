library(dplyr)
library(data.table)
library(metafor)
library(survival)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

mom_list= list()
dad_list= list()
fet_list= list()
q1_list= list()
input= unlist(snakemake@input)

for (coh in cohorts){
input_coh= input[grep(coh, input)]

q1= fread(input_coh[grep('q1', input_coh)])


if (coh== 'harvestm12' | coh== 'harvestm24'){
q1= rename(q1, PREG_ID= PREG_ID_1724)


} else {
q1= rename(q1, PREG_ID= PREG_ID_315)

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

q1= select(q1, PREG_ID, income, edu)

q1$cohort= coh


q1_list[[coh]]= q1

}

fhommom= fread(input[grep('excess_maternal', input)])
fhomdad= fread(input[grep('excess_paternal', input)])
fhomfet= fread(input[grep('excess_fetal', input)])

fhommom= select(fhommom, IID, none_F, cohort)
fhomdad= select(fhomdad, IID, none_F, cohort)
fhomfet= select(fhomfet, IID, none_F, cohort)

names(fhommom)= c('IID', 'FHOM', 'cohort')
names(fhomdad)= c('IID', 'FHOM', 'cohort')
names(fhomfet)= c('IID', 'FHOM', 'cohort')



q1= do.call('rbind', q1_list)

moms= fread(input[grep('mfr_maternal', input)])
dads= fread(input[grep('mfr_paternal', input)])
fets= fread(input[grep('mfr_fetal', input)])


mom= left_join(moms, q1, by= c('PREG_ID', 'cohort')) %>% left_join(., fhommom, by= c('IID', 'cohort'))
dad= left_join(dads, q1, by= c('PREG_ID', 'cohort')) %>% left_join(., fhomdad, by= c('IID', 'cohort'))
fet= left_join(fets, q1, by= c('PREG_ID', 'cohort')) %>% left_join(., fhomfet, by= c('IID', 'cohort'))

mom$FKB= mom$FKB * 100
dad$FKB= dad$FKB * 100
fet$FKB= fet$FKB * 100

mom$KBAVG= mom$KBAVG / 10**6 * 1000
dad$KBAVG= dad$KBAVG / 10**6 * 1000
fet$KBAVG= fet$KBAVG / 10**6 * 1000

mom$tmrca= 100 / (2 * mom$KBAVG)
dad$tmrca= 100 / (2 * dad$KBAVG)
fet$tmrca= 100 / (2 * fet$KBAVG)

mom$tmrca= ifelse(mom$tmrca== Inf, NA, mom$tmrca)
dad$tmrca= ifelse(dad$tmrca== Inf, NA, dad$tmrca)
fet$tmrca= ifelse(fet$tmrca== Inf, NA, fet$tmrca)


exposure= c('FKB', 'KBAVG', 'NSEG', 'FHOM', 'tmrca')


df_l= lapply(setNames(exposure, exposure), function(expo){
surv_formula= as.formula(paste('Surv(SVLEN_UL_DG, spont)~', expo, '+ cohort'))
x= survreg(formula= surv_formula, mom, na.action= na.omit, dist= 'weibull')
return(c(summary(x)$table[2,], n= length(x$linear.predictors), event=sum(!grepl("+", x$y, fixed= T))))})

moms_res= as.data.frame(do.call('rbind', df_l))
moms_res$exposure= rownames(moms_res)
names(moms_res)= c('beta', 'se', 'z', 'pvalue', 'n',  'event','exposure')

moms_res$sample= 'Mothers'
moms_res$model= 'crude'


df_l= lapply(setNames(exposure, exposure), function(expo){
surv_formula= as.formula(paste('Surv(SVLEN_UL_DG, spont)~', expo, '+ cohort'))
x= survreg(formula= surv_formula, dad, na.action= na.omit, dist= 'weibull')
return(c(summary(x)$table[2,], n= length(x$linear.predictors), event=sum(!grepl("+", x$y, fixed= T))))})

dads_res= as.data.frame(do.call('rbind', df_l))
dads_res$exposure= rownames(dads_res)
names(dads_res)= c('beta', 'se', 'z', 'pvalue', 'n',  'event','exposure')


dads_res$sample= 'Fathers'
dads_res$model= 'crude'

df_l= lapply(setNames(exposure, exposure), function(expo){
surv_formula= as.formula(paste('Surv(SVLEN_UL_DG, spont)~', expo, '+ cohort'))
x= survreg(formula= surv_formula, fet, na.action= na.omit, dist= 'weibull')
return(c(summary(x)$table[2,], n= length(x$linear.predictors), event=sum(!grepl("+", x$y, fixed= T))))})

fets_res= as.data.frame(do.call('rbind', df_l))
fets_res$exposure= rownames(fets_res)
names(fets_res)= c('beta', 'se', 'z', 'pvalue', 'n',  'event','exposure')


fets_res$sample= 'Fetuses'
fets_res$model= 'crude'

res_crude= do.call('rbind', list(moms_res, dads_res, fets_res))


df_l= lapply(setNames(exposure, exposure), function(expo){
surv_formula= as.formula(paste('Surv(SVLEN_UL_DG, spont)~', expo, '+ cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +PC10 + income + edu'))
x= survreg(formula= surv_formula, mom, na.action= na.omit, dist= 'weibull')
return(c(summary(x)$table[2,], n= length(x$linear.predictors), event=sum(!grepl("+", x$y, fixed= T))))})

moms_res= as.data.frame(do.call('rbind', df_l))
moms_res$exposure= rownames(moms_res)
names(moms_res)= c('beta', 'se', 'z', 'pvalue', 'n',  'event','exposure')


moms_res$sample= 'Mothers'
moms_res$model= 'adjusted'


df_l= lapply(setNames(exposure, exposure), function(expo){
surv_formula= as.formula(paste('Surv(SVLEN_UL_DG, spont)~', expo, '+ cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +PC10 + income + edu'))
x= survreg(formula= surv_formula, dad, na.action= na.omit, dist= 'weibull')
return(c(summary(x)$table[2,], n= length(x$linear.predictors), event=sum(!grepl("+", x$y, fixed= T))))})

dads_res= as.data.frame(do.call('rbind', df_l))
dads_res$exposure= rownames(dads_res)
names(dads_res)= c('beta', 'se', 'z', 'pvalue', 'n',  'event','exposure')


dads_res$sample= 'Fathers'
dads_res$model= 'adjusted'

df_l= lapply(setNames(exposure, exposure), function(expo){
surv_formula= as.formula(paste('Surv(SVLEN_UL_DG, spont)~', expo, '+ cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +PC10 + income + edu'))
x= survreg(formula= surv_formula, fet, na.action= na.omit, dist= 'weibull')
return(c(summary(x)$table[2,], n= length(x$linear.predictors), event=sum(!grepl("+", x$y, fixed= T))))})

fets_res= as.data.frame(do.call('rbind', df_l))
fets_res$exposure= rownames(fets_res)
names(fets_res)= c('beta', 'se', 'z', 'pvalue', 'n',  'event','exposure')


fets_res$sample= 'Fetuses'
fets_res$model= 'adjusted'

res_adjusted= do.call('rbind', list(moms_res, dads_res, fets_res))

all_res= rbind(res_crude, res_adjusted)

all_res$lowerci = (exp((-1.96* all_res$se) + all_res$beta) -1) * 100
all_res$upperci = (exp((1.96* all_res$se) + all_res$beta) -1) * 100
all_res$HR= (exp(all_res$beta) - 1) * 100

write.table(all_res, snakemake@output[[1]], sep= '\t', row.names=F, col.names= T, quote=F)

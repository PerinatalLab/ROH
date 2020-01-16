library(data.table)
library(dplyr)

df_struct= function(member) {
print(class(unlist(coh_input[grep(paste0('mfr_', member), coh_input)])))
print(unlist(coh_input[grep(paste0('mfr_', member), coh_input)]))
df= fread(unlist(coh_input[grep(paste0('mfr_', member), coh_input)]))


fhom= fread(unlist(coh_input[grep(paste0(member, '_excess'), coh_input)]))
fhom= select(fhom, IID, none_F)
names(fhom)= c('IID', 'FHOM')
df= full_join(df, fhom, by= 'IID')

df$KBAVG= df$KBAVG / 10**6 * 1000
df$tmrca= 100 / (2 * df$KBAVG)
df$tmrca= ifelse(df$tmrca== Inf, NA, df$tmrca)
df$cohort= coh
df$fam= member
df$FKB= df$FKB * 100
df= select(df, FKB, NSEG, FHOM, KBAVG, tmrca, cohort, fam)
return(df)
}

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

mom_list= list()
dad_list= list()
fet_list= list()

input= snakemake@input

for (coh in cohorts){
coh_input= input[grep(coh, input)]

mom= df_struct('maternal')
dad= df_struct('paternal')
fet= df_struct('fetal')

mom_list= c(mom_list, list(mom))
dad_list= c(dad_list, list(dad))
fet_list= c(fet_list, list(fet))
}

mom= do.call('rbind', mom_list)
dad= do.call('rbind', dad_list)
fet= do.call('rbind', fet_list)


sum_mom= mom %>% group_by(cohort) %>% summarize(
FKB_m= median(FKB, na.rm= T), FKB_25= quantile(FKB, probs= 0.25, na.rm=T), FKB_75= quantile(FKB, probs= 0.75, na.rm=T),
NSEG_m= median(NSEG, na.rm= T), NSEG_25= quantile(NSEG, probs= 0.25, na.rm=T), NSEG_75= quantile(NSEG, probs= 0.75, na.rm=T),
FHOM_m= median(FHOM, na.rm= T), FHOM_25= quantile(FHOM, probs= 0.25, na.rm=T), FHOM_75= quantile(FHOM, probs= 0.75, na.rm=T),
KBAVG_m= median(KBAVG, na.rm= T), KBAVG_25= quantile(KBAVG, probs= 0.25, na.rm=T), KBAVG_75= quantile(KBAVG, probs= 0.75, na.rm=T),
tmrca_m= median(tmrca, na.rm= T), tmrca_25= quantile(tmrca, probs= 0.25, na.rm=T), tmrca_75= quantile(tmrca, probs= 0.75, na.rm=T)
)

sum_dad= dad %>% group_by(cohort) %>% summarize(
FKB_m= median(FKB, na.rm= T), FKB_25= quantile(FKB, probs= 0.25, na.rm=T), FKB_75= quantile(FKB, probs= 0.75, na.rm=T),
NSEG_m= median(NSEG, na.rm= T), NSEG_25= quantile(NSEG, probs= 0.25, na.rm=T), NSEG_75= quantile(NSEG, probs= 0.75, na.rm=T),
FHOM_m= median(FHOM, na.rm= T), FHOM_25= quantile(FHOM, probs= 0.25, na.rm=T), FHOM_75= quantile(FHOM, probs= 0.75, na.rm=T),
KBAVG_m= median(KBAVG, na.rm= T), KBAVG_25= quantile(KBAVG, probs= 0.25, na.rm=T), KBAVG_75= quantile(KBAVG, probs= 0.75, na.rm=T),
tmrca_m= median(tmrca, na.rm= T), tmrca_25= quantile(tmrca, probs= 0.25, na.rm=T), tmrca_75= quantile(tmrca, probs= 0.75, na.rm=T)
)

sum_fet= fet %>% group_by(cohort) %>% summarize(
FKB_m= median(FKB, na.rm= T), FKB_25= quantile(FKB, probs= 0.25, na.rm=T), FKB_75= quantile(FKB, probs= 0.25, na.rm=T),
NSEG_m= median(NSEG, na.rm= T), NSEG_25= quantile(NSEG, probs= 0.25, na.rm=T), NSEG_75= quantile(NSEG, probs= 0.25, na.rm=T),
FHOM_m= median(FHOM, na.rm= T), FHOM_25= quantile(FHOM, probs= 0.25, na.rm=T), FHOM_75= quantile(FHOM, probs= 0.25, na.rm=T),
KBAVG_m= median(KBAVG, na.rm= T), KBAVG_25= quantile(KBAVG, probs= 0.25, na.rm=T), KBAVG_75= quantile(KBAVG, probs= 0.25, na.rm=T),
tmrca_m= median(tmrca, na.rm= T), tmrca_25= quantile(tmrca, probs= 0.25, na.rm=T), tmrca_75= quantile(tmrca, probs= 0.25, na.rm=T)
)


d= rbind(mom, dad, fet)

sum_all= d %>% group_by(fam) %>% summarize(
FKB_m= median(FKB, na.rm= T), FKB_25= quantile(FKB, probs= 0.25, na.rm=T), FKB_75= quantile(FKB, probs= 0.75, na.rm=T),
NSEG_m= median(NSEG, na.rm= T), NSEG_25= quantile(NSEG, probs= 0.25, na.rm=T), NSEG_75= quantile(NSEG, probs= 0.75, na.rm=T),
FHOM_m= median(FHOM, na.rm= T), FHOM_25= quantile(FHOM, probs= 0.25, na.rm=T), FHOM_75= quantile(FHOM, probs= 0.75, na.rm=T),
KBAVG_m= median(KBAVG, na.rm= T), KBAVG_25= quantile(KBAVG, probs= 0.25, na.rm=T), KBAVG_75= quantile(KBAVG, probs= 0.75, na.rm=T),
tmrca_m= median(tmrca, na.rm= T), tmrca_25= quantile(tmrca, probs= 0.25, na.rm=T), tmrca_75= quantile(tmrca, probs= 0.75, na.rm=T)
)

write.table(sum_mom, snakemake@output[[1]], col.names= T, row.names= F, sep= '\t', quote=F)
write.table(sum_dad, snakemake@output[[2]], col.names= T, row.names= F, sep= '\t', quote=F)
write.table(sum_fet, snakemake@output[[3]], col.names= T, row.names= F, sep= '\t', quote=F)
write.table(sum_all, snakemake@output[[4]], col.names= T, row.names= F, sep= '\t', quote=F)


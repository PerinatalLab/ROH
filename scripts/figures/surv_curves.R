library(data.table)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)
library(survival)
library(flexsurv)
library(Cairo)

input= unlist(snakemake@input)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')
colors_3= c('#FFBD01', '#00B25D', '#9C02A7')
df_list= list()
dists_long= c('Cohort1', 'Cohort2', 'Cohort3', 'Cohort4', 'Cohort5', 'Cohort6')

parametric_haz <- vector(mode = "list", length = length(cohorts))
n_dists <- length(cohorts)




for (coh in cohorts){
mom= fread(input[grep(coh, input)])

if (coh== 'harvestm12' | coh== 'harvestm24'){
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
names(mom)[names(mom) == 'SentrixID'] <- 'SentrixID_1'
}

mom$cohort= coh
mom= mom[!duplicated(mom$SentrixID_1),]

mom= select(mom, spont, SVLEN_UL_DG, cohort)
df_list= c(df_list, list(mom))

}

d= do.call('rbind', df_list)


for (i in 1:length(cohorts)){
  fit <- flexsurvreg(Surv(SVLEN_UL_DG, spont) ~ 1, data = d[d$cohort== cohorts[i],], dist = 'Weibull') 
  parametric_haz[[i]] <- summary(fit, type = "hazard", ci = FALSE, tidy = TRUE)
  parametric_haz[[i]]$cohort <- dists_long[i]
}

parametric_haz <- rbindlist(parametric_haz)

p1= ggplot(parametric_haz, aes(x = time, y = est, col = cohort)) +
  geom_line(aes(linetype= cohort)) +
theme_cowplot(12, font_size= 12) +
  xlab("Days") + ylab("Hazard") +
  scale_color_manual(name = 'Cohort', values= rep(colors_3, 2)) +
scale_linetype_manual(name= 'Cohort', values= c('solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed')) +
theme(legend.position="bottom")

save_plot(snakemake@output[[1]], plot= p1, device= cairo_ps)

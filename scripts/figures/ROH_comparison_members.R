library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(circlize)
library(viridis)
library(gridExtra)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

df_list= list()


input= unlist(snakemake@input)










for (coh in cohorts){
input_coh= input[grep(coh, input)]

trio= fread(input_coh[grep('trio', input_coh)])


mom= fread(input_coh[grep('mfr_maternal', input_coh)])
dad= fread(input_coh[grep('mfr_paternal', input_coh)])
fet= fread(input_coh[grep('mfr_fetal', input_coh)])

mom= select(mom, NSEG, FKB, KBAVG, IID)
dad= select(dad, NSEG, FKB, KBAVG, IID)
fet= select(fet, NSEG, FKB, KBAVG, IID)

fhommom= fread(input_coh[grep('maternal_excess', input_coh)])
fhomdad= fread(input_coh[grep('paternal_excess', input_coh)])
fhomfet= fread(input_coh[grep('fetal_excess', input_coh)])

fhommom= select(fhommom, IID, none_F)
fhomdad= select(fhomdad, IID, none_F)
fhomfet= select(fhomfet, IID, none_F)

names(fhommom)= paste0(names(fhommom), '_mom')
names(fhomdad)= paste0(names(fhomdad), '_dad')
names(fhomfet)= paste0(names(fhomfet), '_fet')

fhomdf= full_join(trio, fhommom, by= c('Mother'= 'IID_mom')) %>% full_join(., fhomdad, by= c('Father'= 'IID_dad')) %>% full_join(., fhomfet, by= c('Child'= 'IID_fet'))

fhomdf= select(fhomdf, Mother, Father, Child, none_F_mom, none_F_dad, none_F_fet)

names(mom)= paste0(names(mom), '_mom')
names(dad)= paste0(names(dad), '_dad')
names(fet)= paste0(names(fet), '_fet')

df= full_join(trio, mom, by= c('Mother'= 'IID_mom')) %>% full_join(., dad, by= c('Father'= 'IID_dad')) %>% full_join(., fet, by= c('Child'= 'IID_fet'))

df= full_join(df, fhomdf, by= c('Mother', 'Father', 'Child'))

df$cohort= coh

df= select(df, NSEG_mom, FKB_mom, KBAVG_mom, NSEG_dad, FKB_dad, KBAVG_dad, NSEG_fet, FKB_fet, KBAVG_fet, none_F_mom, none_F_dad, none_F_fet, cohort)

#print(summary(df))


df_list= c(df_list, list(df))

}

d= do.call('rbind', df_list)

full_matcor= lapply(cohorts, function(x) { cormat= melt(round(cor(d[d$cohort== x , 1:(ncol(d)-1)], use= 'complete', method= 'spearman'), 2))
cormat$cohort= x
return(cormat)})

matcor= do.call('rbind', full_matcor)

#matcor= melt(cor(d[,1:(ncol(d)-1)], use= 'complete', method= 'spearman'))


level_order= c("FKB_fet"= "Autozygosity offspring", "KBAVG_fet"= "ROH average offspring", "NSEG_fet"= "NSEG offspring", "none_F_fet"= "FHOM offspring", "FKB_mom" ="Autozygosity mother", "KBAVG_mom"= "ROH average mother", "NSEG_mom" = "NSEG mother", "none_F_mom"= "FHOM mother", "FKB_dad"= "Autozygosity father", "KBAVG_dad"= "ROH average father", "NSEG_dad"= "NSEG father", "none_F_dad"= "FHOM father")

var_order= c("FKB_fet", "KBAVG_fet", "NSEG_fet", "none_F_fet", "FKB_mom", "KBAVG_mom", "NSEG_mom" , "none_F_mom", "FKB_dad", "KBAVG_dad", "NSEG_dad", "none_F_dad")

matcor$frisk= ifelse(grepl('NSEG', matcor$Var1), 'NSEG', ifelse(grepl('FKB', matcor$Var1), 'Autozygosity', ifelse(grepl('KBAVG', matcor$Var1), 'ROH average', 'FHOM')))

matcor$Var1= factor(matcor$Var1, levels= var_order)
matcor$Var2= factor(matcor$Var2, levels= var_order)
matcor$cohort= factor(matcor$cohort, levels= cohorts)

matcor$member= ifelse(grepl('mom', matcor$Var1), 'Mother', ifelse(grepl('dad', matcor$Var1), 'Father', 'Offspring'))

s2= ggplot(matcor, aes(frisk, Var2, fill=value))+
  geom_tile(color= "white", size=0.1) + 
facet_grid(cohort ~ member, labeller= labeller(cohort= as_labeller(c('harvestm12'='Cohort1', 'harvestm24'='Cohort2', 'rotterdam1'='Cohort3', 'rotterdam2'='Cohort4', 'normentfeb'='Cohort5', 'normentmay'='Cohort6')))) +
scale_fill_viridis(name = "R coefficient", option ="D") +
scale_x_discrete(labels= level_order) +
scale_y_discrete(labels= level_order) +
theme_bw(base_size=9, base_family = "Source Sans Pro") +
 theme(panel.border = element_blank(),
 panel.grid.major = element_blank(),
 panel.grid.minor = element_blank()) +
 theme(legend.position = "bottom") +
  theme(plot.title=element_text(size = 14)) +
  theme(axis.text.y=element_text(size=9)) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(axis.text=element_text(size=9)) +
#  theme(axis.text.x= element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_text(size=9)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.key.width = unit(1,"cm")) +
  theme(strip.text= element_text(size = 9)) +
  theme(axis.title.x=element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))+
  theme(panel.spacing.x = unit(0, "lines")) +
theme(panel.spacing.x=unit(0.1, "lines"), panel.spacing.y=unit(0.1, "lines"))




d$id= 1:nrow(d)

FKB= d %>% gather(member, value, c('FKB_mom', 'FKB_dad', 'FKB_fet')) %>% select(member, value, id)
FHOM= d %>% gather(member, value, c('none_F_mom', 'none_F_dad', 'none_F_fet')) %>% select(member, value, id)
KBAVG= d %>% gather(member, value, c('KBAVG_mom', 'KBAVG_dad', 'KBAVG_fet')) %>% select(member, value, id)
NSEG= d %>% gather(member, value, c('NSEG_mom', 'NSEG_dad', 'NSEG_fet')) %>% select(member, value, id)

FKB$value= FKB$value * 100
KBAVG$value= KBAVG$value / (10**4 * 100)

x= do.call('bind_rows', list(FKB, FHOM, KBAVG, NSEG))

x$frisk= ifelse(grepl('NSEG', x$member), 'NSEG', ifelse(grepl('FKB', x$member), 'Autozygosity', ifelse(grepl('KBAVG', x$member), 'ROH average', 'FHOM')))

x$member= ifelse(grepl('dad', x$member), 'Father', ifelse(grepl('mom', x$member), 'Mother', 'Offspring'))
x$member <- factor(x$member, levels = c('Mother', 'Offspring', 'Father'))

x$frisk= factor(x$frisk, levels= c('Autozygosity', 'NSEG', 'ROH average', 'FHOM'))

p= ggplot(x, aes(sample = value, colour = factor(member))) +
  stat_qq(size= 1) +
  stat_qq_line() +
scale_colour_viridis_d(name = "Family member", labels = c("Mother", "Offspring", "Father")) +
facet_wrap(vars(frisk), nrow= 2, ncol= 2, scales= 'free', labeller= labeller(frisk= as_labeller(c('Autozygosity'= 'Autozygosity', 'FHOM'='FHOM', 'NSEG'='Number of segments', 'ROH average'='Average segment length')))) +
theme_bw(base_size=9, base_family = "Source Sans Pro") +
 theme(panel.border = element_blank(), 
 panel.grid.major = element_blank(), 
 panel.grid.minor = element_blank()) +
 theme(legend.position = "bottom") +
  theme(plot.title=element_text(size = 14)) +
  theme(axis.text=element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(axis.text=element_text(size=9)) +
  xlab('Theoretical quantiles') +
  ylab('Observed') +
  theme(legend.title=element_text(size=9)) +
  theme(legend.text=element_text(size=9)) +
  theme(strip.text= element_text(size = 9)) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))


grid.arrange
ggsave(snakemake@output[[1]], s2, dpi= 'retina')
ggsave(snakemake@output[[2]], p, dpi= 'retina')

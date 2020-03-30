library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

df_struct= function(member) {
df= fread(unlist(input[grep(member, input)]))
df$KBAVG= df$KBAVG / 10**6 * 1000
df$tmrca= 100 / (2 * df$KBAVG)
df$fam= member
df= select(df, tmrca, cohort, fam)
return(df)
}

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')
colors_3= c('#FFBD01', '#00B25D', '#9C02A7')

mom_list= list()
dad_list= list()
fet_list= list()

input= snakemake@input


mom= df_struct('maternal')
dad= df_struct('paternal')
fet= df_struct('fetal')

d= rbind(mom, dad, fet)

d$cohort= factor(d$cohort, levels= cohorts)
d$fam= ifelse(d$fam== 'maternal', 'Maternal', ifelse(d$fam== 'fetal', 'Fetal', ifelse(d$fam== 'paternal', 'Paternal', NA)))
d$fam= factor(d$fam, levels= c('Maternal', 'Paternal', 'Fetal'))
d$tmrca= ifelse(d$tmrca== Inf, NA , d$tmrca)

p1= ggplot(d, aes(x= fam, y= tmrca, colour= fam, colour= fam)) +
geom_dotplot(binaxis='y', stackdir='center', dotsize= 0.3, binwidth= 0.5) +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(labels= c('maternal'='Maternal', 'paternal'= 'Paternal', 'fetal'= 'Fetal'), values= colors_3) +
 theme(legend.position = "none") +
        xlab('Member') +
        ylab('Time to most recent common ancestor, generations') +
  theme(strip.background = element_rect(colour="white", fill="white"))

save_plot(snakemake@output[[1]], p1, base_width=297, base_height=210, units="mm", device= cairo_ps)



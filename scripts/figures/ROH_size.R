library(data.table)
library(ggridges)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)
library('ggrepel')

colors_6= c('#FFBD01', '#00B25D', '#9C02A7', 'black', 'grey', 'cyan')

d= fread(snakemake@input[[1]])
d$fam= 'Maternal'
d1= fread(snakemake@input[[2]])
d1$fam= 'Paternal'
d2= fread(snakemake@input[[3]])
d2$fam= 'Fetal'

#df= group_by(x, ROH_class) %>% summarize(minROH= min(min_distance, na.rm=T)/ 1000, maxROH= max(min_distance, na.rm=T)/ 1000)

d= do.call('rbind', list(d, d1, d2))

d= mutate(d, cohort= ifelse(cohort =='harvestm12', 'Cohort1', ifelse(cohort== 'harvestm24', 'Cohort2', ifelse(cohort== 'rotterdam1', 'Cohort3', ifelse(cohort== 'rotterdam2', 'Cohort4', ifelse(cohort== 'normentfeb', 'Cohort5', 'Cohort6'))))))


p1= ggplot() + 
#geom_density(data= d, aes(KB/1000, group= cohort, color= cohort)) +
geom_density_ridges(data= d, aes(x= KB/1000, y= cohort, group= cohort, fill= cohort), alpha= 0.5) +
theme_cowplot(12, font_size= 10) +
scale_fill_manual(name = "Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= colors_6) +
facet_wrap(vars(fam), ncol= 3) +
theme(legend.position="none",
#strip.text = element_blank(),
          strip.background = element_blank()) +
#geom_rect(aes(xmin= unlist(df[df$ROH_class== 2, 'minROH'])/1000, xmax = unlist(df[df$ROH_class== 2, 'maxROH'])/1000, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 1/10) +
#geom_rect(aes(xmin= unlist(df[df$ROH_class== 3, 'minROH'])/1000, xmax = unlist(df[df$ROH_class== 3, 'maxROH'])/1000, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 1/10) +
xlab('Autozygous segment size, cM') +
ylab('Density') +
scale_x_continuous(limits= c(0, 10), expand = c(0, 0)) +
theme(panel.spacing = unit(1, "lines"))


save_plot(file= snakemake@output[[1]], plot= p1, base_width= 190, base_height=120, units="mm" , device= cairo_ps)

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')
colors_4= c('#FFBD01', '#00B25D', '#9C02A7', '#000000')

df_list= list()

input= snakemake@input

for (coh in cohorts){
arg= fread(unlist(input[grep(coh, input)]))

arg$fname= lapply(arg$file, function(x) unlist(strsplit(x, '/'))[8])
arg$pruning= ifelse(grepl('hard', arg$fname), 3, ifelse(grepl('moderate', arg$fname), 2, ifelse(grepl('soft', arg$fname), 1,0)))
arg$bp_cm= ifelse(grepl('_bp', arg$fname), 'BP', 'cM')

arg$fname= gsub('.*fetal_', '',arg$fname)
arg$fname= gsub('1e-07', '0.0000001', arg$fname)
arg$fname= gsub('.hom.indiv', '',arg$fname)
arg= arg %>% separate(fname, c('dens', 'SNP', 'length', 'het', 'GAP'), sep= '_')

arg$cohort= coh

df_list= c(df_list, list(arg))
}

d= do.call('rbind', df_list)

sorted_labels <- sort(as.integer(levels(as.factor(d$SNP))))
d$SNP= factor(d$SNP, levels = sorted_labels)
d$cohort= factor(d$cohort, levels= cohorts)

x= d %>% group_by(SNP, bp_cm, pruning, het) %>% summarize(R2= mean(R2))


##### Figure 1.

bp= filter(d, bp_cm== 'BP')
cm= filter(d, bp_cm== 'cM')

colnames(bp)[1] <- "R2_bp"
colnames(cm)[1] <- "R2_cm"
df= inner_join(bp, cm, by= c('SNP', 'het', 'cohort', 'pruning'))

p1= ggplot(x, aes(x= as.numeric(as.character(SNP)), y= R2, colour= as.factor(pruning))) +
geom_point(size= 1) +
geom_line() + 
scale_x_continuous(breaks=c(0, 400, 25), expand = c(0,0)) +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name= expression(paste('Pruning ', R^2, ' threshold')), labels= c('None', '0.9', '0.5', '0.1'), values= colors_4) +
facet_grid(het ~ bp_cm, labeller = labeller(het= as_labeller(c('0'= 'No het. allowed', '1'= 'One het. allowed')),
bp_cm= as_labeller(c('BP'='Physical distance', 'cM'='Genetic distance')))) +
scale_y_continuous(breaks=seq(round(min(x$R2, na.rm= T), 1), 0.2, 0.6), expand = c(0,0)) +
 theme(legend.position = "bottom") +
        xlab('Number of genetic variants included') +
        ylab(expression(paste('Correlation coefficient between offspring ROH and parental genetic relatedness'))) +
  theme(strip.background = element_rect(colour="white", fill="white"))

p2= ggplot(df, aes(x= R2_bp, y= R2_cm, colour= cohort, shape= cohort)) +
geom_point(size= 2) +
scale_x_continuous(expand = c(0,0)) +
scale_y_continuous(expand = c(0,0)) +
geom_rug(col="grey", alpha=0.1, size=1.5)+
geom_abline(intercept = 0, slope = 1) +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name="Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= rep(colors_3,2)) +
scale_shape_manual(name="Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6'), values= rep(15:16, 3)) +
 theme(strip.text= element_text()) +
 xlab(expression(paste('Correlation coefficient between offspring ROH and parental genetic relatedness using physical distance'))) +
        ylab(expression(paste('Correlation coefficient between offspring ROH and parental genetic relatedness using genetic distance')))






save_plot(snakemake@output[[1]], p1, base_width=297, base_height=210, units="mm", device= cairo_ps)
save_plot(snakemake@output[[2]], p2, base_width=297, base_height=210, units="mm", device= cairo_ps)




#### Supp figure



s1B= ggplot(filter(d, bp_cm== 'cM', het== 0), aes(x= as.numeric(as.character(SNP)), y= R2, colour= as.factor(pruning))) +
geom_point(size= 1) +
geom_line() +
scale_x_continuous(breaks= c(0, 400, 25), expand = c(0,0)) +
scale_y_continuous(breaks= seq(round(min(x$R2, na.rm= T), 1), 0.6, 0.2), expand = c(0,0)) +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name= expression(paste('Pruning ', R^2, ' threshold')), labels= c('None', '0.9', '0.5'), values= colors_4) +
facet_wrap(~cohort, labeller = labeller(cohort= as_labeller(c('harvestm12'='Cohort1', 'harvestm24'='Cohort2', 'rotterdam1'='Cohort3', 'rotterdam2'='Cohort4', 'normentfeb'='Cohort5', 'normentmay'='Cohort6'))), nrow=2, ncol=3) +
 theme(legend.position = "bottom") +
        xlab('Number of genetic variants included') +
        ylab(expression(paste('Correlation coefficient between offspring ROH and parental genetic relatedness'))) +
  theme(strip.background = element_rect(colour="white", fill="white"))

s1A= ggplot(filter(d, bp_cm== 'cM', het== 1), aes(x= as.numeric(as.character(SNP)), y= R2, colour= as.factor(pruning))) +
geom_point(size= 1) +
geom_line() +
scale_x_continuous(breaks=c(0, 400, 25), expand = c(0,0)) +
scale_y_continuous(breaks=seq(round(min(x$R2, na.rm= T), 1), 0.6, 0.2), expand = c(0,0)) +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name= expression(paste('Pruning ', R^2, ' threshold')), labels= c('None', '0.9', '0.5'), values= colors_4) +
facet_wrap(~cohort, labeller = labeller(cohort= as_labeller(c('harvestm12'='Cohort1', 'harvestm24'='Cohort2', 'rotterdam1'='Cohort3', 'rotterdam2'='Cohort4', 'normentfeb'='Cohort5', 'normentmay'='Cohort6'))), ncol= 3, nrow= 2) +
 theme(legend.position = "bottom") +
        xlab('Number of genetic variants included') +
        ylab(expression(paste('Correlation coefficient between offspring ROH and parental genetic relatedness'))) +
  theme(strip.background = element_rect(colour="white", fill="white"))


save_plot(snakemake@output[[3]], plot= s1A, base_width=297, base_height=210, units="mm", device= cairo_ps)
save_plot(snakemake@output[[4]], plot= s1B, base_width=297, base_height=210, units="mm", device= cairo_ps)

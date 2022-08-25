library(data.table)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)
library('ggrepel')

colors_3= c('#9C02A7', '#FFBD01', '#00B25D')

d= fread(snakemake@input[[1]], col.names= c('gene','n', 'freq', 'beta', 'se', 'pvalue', 'loglik'))

d= separate(d, gene, into= c('chr', 'gene', 'EntrezID'), sep= ':')
names(d)= paste0(names(d), '_mom')

x= sum(as.numeric(readLines(snakemake@input[[2]])))

df= fread(snakemake@input[[3]], col.names= c('gene','n', 'freq', 'beta', 'se', 'pvalue', 'loglik'))

df= separate(df, gene, into= c('chr', 'gene', 'EntrezID'), sep= ':')

x1= sum(as.numeric(readLines(snakemake@input[[4]])))

d= inner_join(d, df, by= c('gene_mom'= 'gene'))

d$HC= ifelse(d$pvalue< 0.05/x1, 'Fetal', NA)
d$HC= ifelse(d$pvalue_mom< 0.05/x, 'Maternal', NA)
d$HC= ifelse(d$HC== NA, 'None', d$HC)

p1= ggplot(filter(d, pvalue< 0.05 | pvalue_mom< 0.05)) +
   geom_point(aes(x=beta_mom, y= beta, colour= HC, fill= HC), size=0.3, shape= 21) +   # Show all points
theme_cowplot(12, font_size= 12) + #theme_minimal_hgrid(12, rel_small = -1) + 
scale_colour_manual(values= c('Maternal'= colors_3[1], 'Fetal'= colors_3[2], 'None'= 'grey'), labels= c('Maternal', 'Fetal', 'None')) +
scale_fill_manual(values= c('Maternal'= colors_3[1], 'Fetal'= colors_3[2], 'None'= 'grey'), labels= c('Maternal', 'Fetal', 'None')) +
xlab('Maternal effect size') +
ylab('Fetal effect size') +
geom_hline(yintercept= 0, size= 0.2, colour= 'black', alpha= 0.01) +
geom_vline(xintercept= 0, size= 0.2, colour= 'black', alpha= 0.01)

save_plot(file= snakemake@output[[1]], plot= p1, base_width=120, base_height=120, units="mm", device= cairo_ps)



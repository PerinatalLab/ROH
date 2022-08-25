library('data.table')
library('ggplot2')
library('dplyr')
library('tidyr')
library('ggrepel')
library('cowplot')
library('survminer')
library(flexsurv)

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')

d= fread(snakemake@input[[1]])
d$fam= 'Maternal'

d1= fread(snakemake@input[[2]])
d1$fam= 'Paternal'

d2= fread(snakemake@input[[3]])
d2$fam= 'Fetal'

d= do.call('rbind', list(d, d1, d2))


x1= ggplot(d, aes(rel_overlap, group= fam, colour= fam)) + 
geom_density() +
theme(legend.position = 'bottom') +
theme_cowplot(12, font_size= 12) +
scale_colour_manual(name = "Sample", values= colors_3) +
xlab('Relative segment overlap') +
    ylab('Density')


save_plot(snakemake@output[[1]], plot= x1, base_width=120, base_height=100, units="mm")


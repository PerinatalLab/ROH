library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)

colors_3= c('#9C02A7', '#FFBD01', '#00B25D')


input= snakemake@input

d= fread(unlist(input[grep('maternal', input)]))
d$fam= 'Maternal'

d1= fread(unlist(input[grep('paternal', input)]))
d1$fam= 'Paternal'

d2= fread(unlist(input[grep('fetal', input)]))
d2$fam= 'Fetal'

d= do.call('rbind', list(d, d1, d2))

d= filter(d, !is.na(pvalue))
df= arrange(d, pvalue) %>% group_by(fam) %>% mutate(exp1= -log10(1:length (pvalue)/length (pvalue)))

p1= ggplot(df, aes(exp1, -log10(pvalue), colour= factor(fam))) + 
  geom_point(size= 0.4) + 
  geom_abline(intercept = 0, slope = 1, alpha = .5) +
scale_colour_manual(values= colors_3, name='') +
theme_cowplot(12, font_size= 12) +
xlab('Expected (-log10(p-value))') +
ylab('Observed (-log10(p-value))')


save_plot(file= snakemake@output[[1]], plot= p1, base_width=120, base_height=120, units="mm", device= cairo_ps)

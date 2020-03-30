  library('data.table')
  library('ggplot2')
  library('dplyr')
  library('tidyr')
  library('ggrepel')
  library('cowplot')
  

d= fread(snakemake@input[[1]], h=T)
geno= fread(snakemake@input[[2]], h=T)
cm= fread(snakemake@input[[3]], h=T)


names(cm)= c('chr', 'pos', 'rate', 'cM')
geno$start= geno$cM1 * 10**6
geno$end= geno$cM2 * 10**6


if (nrow(d)==0){
p1= ggplot() + theme_void()
save_plot(snakemake@output[[1]], plot= p1, base_width=297, base_height=210, units="mm", device= cairo_ps)

} else {
#chrs= unique(d[1, 'chr'])

df= fread(snakemake@input[[4]])
df= mutate(df, chr= as.numeric(chr), cM1= as.numeric(cM1), cM2= as.numeric(cM2))

df= arrange(df, pvalue)
chrs= unique(df[1, 'chr'])

df= filter(df, !duplicated(segment), chr %in% chrs, !is.na(pvalue))

df$labs= paste('Chromosome', df$chr)

cm_range= 1

cM_df= group_by(df, chr) %>% do(head(., 1)) %>% summarize(cMm= (cM1 + cM2) / 2)
num_cols= ifelse(length(chrs)<= 3, length(chrs), round(length(chrs) / 2))

df= inner_join(cM_df, df, by= 'chr')

df$cM1= ifelse((df$cM1 < df$cMm - cm_range*10**6) & (df$cM2<= df$cMm + cm_range*10**6), df$cMm - cm_range*10**6, df$cM1)
df$cM2= ifelse((df$cM1>= df$cMm - cm_range*10**6) & (df$cM2> df$cMm + cm_range*10**6), df$cMm + cm_range*10**6, df$cM2)

df= df %>% filter(cM1>= cMm - cm_range*10**6, cM2<= cMm + cm_range*10**6)

x= group_by(df, chr) %>% summarize(cM1= min(cM1), cM2= max(cM2))
x= arrange(x, chr, cM1)
x= group_by(x, chr) %>% mutate(empt= cM1 - shift(cM2))
x$start= ifelse(x$empt> 0, shift(x$cM2), NA)
x$end= ifelse(x$empt> 0, x$cM1, NA)
#x= filter(x, !is.na(start))

if (nrow(x)>0 ){
x= group_by(x, chr) %>% filter(!duplicated(start))

x$labs= paste('Chromosome', x$chr)
}

newdf= group_by(df, chr) %>% summarize(cMm= unique(cMm))
newdf= mutate(newdf, cM1= cMm  - (cm_range * 10**6), cM2= cMm + (cm_range * 10**6))
newdf$labs= paste('Chromosome', newdf$chr)

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')
colors_2= c('#00B25D', '#9C02A7')

cM1= df[df$pvalue== min(df$pvalue), 'cM1']
cM2= df[df$pvalue== min(df$pvalue), 'cM2']

cm= filter(cm, chr %in% chrs)
cm= inner_join(cm, newdf[, c('chr', 'cM1', 'cM2')], by= 'chr')
cm= filter(cm, cM >= cM1 /10**6, cM<= cM2/10**6)



geno= group_by(geno, chr, gene) %>% summarize(start= min(start), end= max(end)) %>% filter(!duplicated(gene))
geno$labs= paste('Chromosome', geno$chr)
geno= inner_join(geno, newdf[, c('chr', 'cM1', 'cM2')], by= 'chr')
geno= geno %>% rowwise() %>% mutate(overlap= max(0, min(end, cM2) - max(start, cM1))) %>% filter(overlap> 0)

geno$rk= with(geno, rank(end - start))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

CHR== chrs


dff= filter(df, chr== CHR)
newdff= filter(newdf, chr== CHR)
xf= filter(x, chr== CHR)
cmf= filter(cm, chr==CHR)
genof= filter(geno, chr== CHR)
genof= arrange(genof, desc(rk))
genof= genof[1:10, ]
genof= filter(genof, !is.na(chr))
p1= ggplot() +
geom_step(data=cmf, aes(cM, rate ), colour= 'black', alpha= 0.7) +
    theme_cowplot(12) +
 geom_segment(data= newdff, aes(x= cM1/10**6, xend= cM2/10**6, y= 1, yend= 1), alpha= 0) +
    facet_wrap(~ labs, scales = "free_x", ncol= num_cols) +
    scale_colour_manual(values= colors_2) +
    theme(strip.background = element_blank(),
	strip.text = element_blank(),
          legend.position="none") +
    ylab('Recombination rate, cM/bp') +
    xlab('Genetic distance, cM') + 
ylim(0, 100)+
coord_cartesian(expand = FALSE, xlim= c(min(newdff$cM1) /10**6, max(newdff$cM2) /10**6))

p2= ggplot() +
 geom_segment(data= dff, aes(x=cM1 / 10**6, xend=cM2/ 10**6, y= -log10(pvalue), yend= -log10(pvalue),colour= as.factor(as.numeric(beta>0)))) +
    scale_colour_manual(values= colors_2) +
    theme_cowplot(12) +
 geom_segment(data= newdff, aes(x= cM1/10**6, xend= cM2/10**6, y= 1, yend= 1), alpha= 0) +
    facet_wrap(~ labs, scales = "free_x", ncol= num_cols) +
    geom_segment(data= xf, aes(x= start/10**6, xend= end/10**6, y= 0, yend=0), colour= 'black') +
    theme(strip.background = element_blank(),
          legend.position="none",
	strip.text = element_text(size= 12),
	axis.line.x=element_blank(),
	axis.title.x=element_blank(),
	axis.text.x=element_blank(),
	axis.ticks.x=element_blank()) +
    ylab('Pvalue, -log10') +
    xlab('Genetic distance, cM') +
coord_cartesian(ylim= c(0, max(-log10(dff$pvalue)) + 1) , expand = FALSE, xlim= c(min(newdff$cM1) /10**6, max(newdff$cM2) /10**6))

p3= ggplot() +
   facet_wrap(~ labs, scales = "free_x", ncol= num_cols) +
 geom_segment(data= genof, aes(x= start/10**6, xend= end/10**6, y= range01(rk), yend= range01(rk)), colour= '#9C02A7') +
  geom_text_repel(data= genof, aes(x= (start/10**6 + end/10**6)/2, y= range01(rk), label= gene), size= 3,  hjust = 0.5, force= 1, vjust= 1, colour= 'black') +
geom_segment(data= newdff, aes(x= cM1/10**6, xend= cM2/10**6, y= 1, yend= 1), alpha= 0) +
    facet_wrap(~ labs, scales = "free_x", ncol= num_cols) +
    theme_cowplot(12) +
    coord_cartesian( ylim = c(0, 1.01), expand = FALSE, xlim= c(min(newdff$cM1) /10**6, max(newdff$cM2) /10**6)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())

x_plot= plot_grid(p2, p1, p3, align= 'v',  ncol= 1, rel_heights= c(5,4,2))




save_plot(snakemake@output[[1]], plot= x_plot, base_width=297, base_height=210, units="mm", device= cairo_ps)

}

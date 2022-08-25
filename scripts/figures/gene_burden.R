library('data.table')
  library('ggplot2')
  library('dplyr')
  library('tidyr')
  library('ggrepel')
  library('cowplot')








geno= fread(snakemake@input[[1]], h=T)
cm= fread(snakemake@input[[2]], h=T)
cd= fread(snakemake@input[[4]], h=T)

names(cm)= c('chr', 'pos', 'rate', 'cM')


cd= cd[, c('chr', 'EntrezID', 'cdcM1', 'cdcM2', 'strand')]
cd$cdcM1= cd$cdcM1 * 10**6
cd$cdcM2= cd$cdcM2 * 10**6

geno= inner_join(geno, cd, by= c('EntrezID', 'chr'))

txstart= ifelse(geno$strand== '-', geno$cdcM2, geno$cdcM1)
txend= ifelse(geno$strand== '-', geno$cdcM1, geno$cdcM2)

geno$cdcM1= txstart
geno$cdcM2= txend

geno$start= geno$cM1 * 10**6
geno$end= geno$cM2 * 10**6

df= fread(snakemake@input[[3]])
names(df)= c('symbol', 'samplesize', 'freq', 'beta', 'se', 'pvalue', 'loglik')
df= separate(df, symbol, into= c('chr', 'gene', 'EntrezID'), sep= ':')
df= mutate(df, chr= as.numeric(chr), EntrezID= as.numeric(EntrezID))

df= arrange(df, pvalue)
chrs= unique(df[1, 'chr'])


df= filter(df, !duplicated(EntrezID), chr %in% chrs, !is.na(pvalue))

df$labs= paste('Chromosome', df$chr)

cm_range= 1

df= inner_join(df, geno, on= 'EntrezID')

cM_df= group_by(df, chr) %>% do(head(., 1)) %>% summarize(cMm= (start + end) / 2)
num_cols= ifelse(length(chrs)<= 3, length(chrs), round(length(chrs) / 2))

df= inner_join(cM_df, df, by= 'chr')

df$start= ifelse((df$start < df$cMm - cm_range*10**6) & (df$end<= df$cMm + cm_range*10**6), df$cMm - cm_range*10**6, df$start)
df$end= ifelse((df$start>= df$cMm - cm_range*10**6) & (df$end> df$cMm + cm_range*10**6), df$cMm + cm_range*10**6, df$end)

df= df %>% filter(start>= cMm - cm_range*10**6, end<= cMm + cm_range*10**6)

x= group_by(df, chr) %>% summarize(cM1= min(start)* 10**6, cM2= max(end)* 10**6)
x= arrange(x, chr, cM1)
x= group_by(x, chr) %>% mutate(empt= cM1 - shift(cM2))
x$start= ifelse(x$empt> 0, shift(x$cM2), NA)
x$end= ifelse(x$empt> 0, x$cM1, NA)

x$labs= paste('Chromosome', x$chr)


newdf= group_by(df, chr) %>% summarize(cMm= unique(cMm))
newdf= mutate(newdf, cM1= cMm  - (cm_range * 10**6), cM2= cMm + (cm_range * 10**6))
newdf$labs= paste('Chromosome', newdf$chr)

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')
colors_2= c('#00B25D', '#9C02A7')

cM1= df[df$pvalue== min(df$pvalue), 'start']
cM2= df[df$pvalue== min(df$pvalue), 'end']

cm= filter(cm, chr %in% chrs)
cm= inner_join(cm, newdf[, c('chr', 'cM1', 'cM2')], by= 'chr')
cm= filter(cm, cM >= cM1 /10**6, cM<= cM2/10**6)


geno= group_by(geno, chr, gene) %>% summarize(start= min(start), end= max(end), cdcM1= min(cdcM1), cdcM2= max(cdcM2)) %>% filter(!duplicated(gene))
geno$labs= paste('Chromosome', geno$chr)
geno= inner_join(geno, newdf[, c('chr', 'cM1', 'cM2')], by= 'chr')
geno= geno %>% rowwise() %>% mutate(overlap= max(0, min(end, cM2) - max(start, cM1))) %>% filter(overlap> 0)

geno$rk= with(geno, rank(end - start))

range01 <- function(x){
if (length(x)>1){
(x-min(x))/(max(x)-min(x))
} else{
1
}
}

CHR= chrs


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
 geom_segment(data= dff, aes(x=start / 10**6, xend=end/ 10**6, y= -log10(pvalue), yend= -log10(pvalue),colour= as.factor(as.numeric(beta>0)))) +
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

p3=ggplot() +
facet_wrap(~ labs, scales = "free_x", ncol= num_cols) +
 geom_segment(data= genof, aes(x= start/10**6, xend= end/10**6, y= range01(rk), yend= range01(rk)), colour= 'black') +
 geom_segment(data= genof, aes(x= cdcM1/10**6, xend= cdcM2/10**6, y= range01(rk), yend= range01(rk)), colour= '#9C02A7', size= 1.5, alpha= 0.9, arrow = arrow(length = unit(0.1, "cm")), lineend = "butt", linejoin = "mitre") +
  geom_text_repel(data= genof, aes(x= end/10**6, y= range01(rk), label= gene), size= 3,  hjust = 1, force= 1, vjust= 0.5, colour= 'black') +
geom_segment(data= newdff, aes(x= cM1/10**6, xend= cM2/10**6, y= 1, yend= 1), alpha= 0) +
    facet_wrap(~ labs, scales = "free_x", ncol= num_cols) +
    theme_cowplot(12) +
    coord_cartesian( ylim = c(-0.1, 1.1), expand = FALSE, xlim= c(min(newdff$cM1) /10**6, max(newdff$cM2) /10**6)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())

x_plot= plot_grid(p2, p1, p3, align= 'v',  ncol= 1, rel_heights= c(5,4,2))




save_plot(snakemake@output[[1]], plot= x_plot, base_width=297, base_height=210, units="mm", device= cairo_ps)


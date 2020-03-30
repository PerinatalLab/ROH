library(data.table)
library(ggrepel)
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

df= fread(snakemake@input[[1]])

chrs= unique(df$chr)

geno= fread(snakemake@input[[3]], h=T)
geno$start= geno$cM1 * 10**6
geno$end= geno$cM2 * 10**6

geno= group_by(geno, chr, gene) %>% summarize(start= min(start), end= max(end)) %>% filter(!duplicated(gene))

if (nrow(df)==0){
p1= ggplot() + theme_void()
save_plot(snakemake@output[[1]], plot= p1, base_width=297, base_height=210, units="mm", device= cairo_ps)
} else {

df= separate(df, segment, into= c ('chr', 'cM1', 'cM2'), sep=':')

df= mutate(df, chr= as.numeric(chr), cM1= as.numeric(cM1), cM2= as.numeric(cM2))

df_list= list()

for (coh in cohorts){
coh_input= input[grep(coh, input)]

hom= fread(unlist(coh_input[grepl('hom', coh_input)]))


hom= inner_join(hom, df, by= c('CHR' = 'chr'))

hom$labs= paste('Chromosome', hom$CHR)
hom$cohort= coh
hom= hom %>% rowwise() %>% mutate(overlap= max(0, min(POS2, cM2) - max(POS1, cM1)))
hom= hom[!duplicated(hom[, c('IID', 'POS1', 'POS2')]),]

hom= filter(hom, overlap> 0) %>% select( POS1, POS2, labs, IID, cohort, CHR)

df_list[[coh]]= hom

}

d= do.call('rbind', df_list)


pheno= fread(snakemake@input[[2]])

pheno= select(pheno, IID, SVLEN_UL_DG, spont, cohort)
pheno= filter(pheno, !is.na(SVLEN_UL_DG), !is.na(spont))
d= inner_join(d, pheno, by= c('IID', 'cohort'))

ncolors= length(unique(d$cohort))

if (ncolors<= 3){
color_values=  colors_3
} else if (ncolors==4) {
color_values= c(colors_3, 'black')
} else if (ncolors==5) {
color_values=c(colors_3, 'black', 'grey')
} else {
color_values= c(colors_3, 'black', 'grey', 'cyan')
}

d$survobj= paste0('(', d$SVLEN_UL_DG, '; ', d$spont, ')')
d$IID= paste(d$survobj, d$IID, sep= '_')

newdf= group_by(d, CHR) %>% summarize(cM1min= min(POS1, na.rm=T), cM2max= max(POS2, na.rm=T))

geno$labs= paste('Chromosome', geno$chr)
geno= inner_join(geno, newdf, by= c('chr'= 'CHR'))
geno= mutate(geno, start= ifelse((start < cM1min) & (end<= cM2max), cM1min, start), end= ifelse((start>= cM1min) & (end> cM2max), cM2max, end))

geno= filter(geno, start>= cM1min, end <= cM2max)

#geno= geno %>% rowwise() %>% mutate(overlap= max(0, min(end, cM2) - max(start, cM1))) %>% filter(overlap> 0)

geno$rk= with(geno, rank(end - start))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

plot_individuals= function(chrom){
df= filter(d, CHR== chrom)
genof= filter(geno, chr== chrom)
genof= arrange(genof, desc(rk))
genof= genof[1:10, ]

p1= ggplot() +
  geom_segment(data= df, aes(x= POS1 / 10**6, xend= POS2/ 10**6, y= as.factor(IID), yend= as.factor(IID), colour= cohort), size=1) +
  theme_cowplot(12, font_size= 12) +
facet_wrap(~labs) +
  scale_y_discrete(breaks = d$IID, labels= d$survobj) +
scale_colour_manual(values= color_values) +
theme(strip.background = element_blank(),
	strip.text = element_text(size = 12),
        legend.position="none") +
  ylab('Time-to-spontaneous delivery; Event)') +
  xlab('Genetic distance, cM')

p2= ggplot() +
   facet_wrap(~ labs, scales = "free_x") +
 geom_segment(data= genof, aes(x= start/10**6, xend= end/10**6, y= range01(rk), yend= range01(rk)), colour= '#9C02A7') +
  geom_text_repel(data= genof, aes(x= (start/10**6 + end/10**6)/2, y= range01(rk), label= gene), size= 3,  hjust = 0.5, force= 1, vjust= 1, colour= 'black') +
    theme_cowplot(12, font_size= 12) +
    coord_cartesian( ylim = c(0, 1.01), expand = FALSE, xlim= c(min(genof$cM1min) /10**6, max(genof$cM2max) /10**6)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())

plot_grid(p1, p2, align= 'v', ncol= 1, rel_heights= c(6,3))
}

plots= lapply(chrs, plot_individuals)

if (length(chrs)>1){

x_plot= plot_grid(plots[[1]], plots[[2]], ncol = 1)
} else{
x_plot= plots[[1]]
}

save_plot(snakemake@output[[1]], plot= x_plot, base_width=297, base_height=210, units="mm", device= cairo_ps)

}

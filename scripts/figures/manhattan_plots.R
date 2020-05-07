library(data.table)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)
library('ggrepel')

d= fread(snakemake@input[[1]])
d= separate(d, segment, into= c('chr', 'cM1', 'cM2'), sep= ':')
d$cM1= as.numeric(d$cM1)
d$cM2= as.numeric(d$cM2)
d$chr= as.numeric(d$chr)
d$z= d$beta/ d$sd
df= d

df$mcM= (df$cM1 + df$cM2) / 2
eff= sum(as.numeric(readLines(snakemake@input[[2]])))
#colors_3= c('#FFBD01', '#00B25D', '#9C02A7')
colors_3= c('#9C02A7', '#FFBD01', '#00B25D')

gene= fread(snakemake@input[[3]])


df= df %>% filter(!is.na(chr))
  don <- df %>%
    group_by(chr)      %>%
    summarise(chr_len= max(mcM)) %>%
    mutate(tot= cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(df, ., by= 'chr') %>%
    arrange(chr, mcM) %>% # Add a cumulative position of each SNP
    mutate( BPcum=mcM+tot)

  axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('chr', 'center')

HC= 0.05/eff
gene= inner_join(select(don, segment, BPcum), gene, by= 'segment')
moms= ggplot(don) +
    geom_point(aes(x=BPcum, y= z, colour= -log10(pvalue), fill= -log10(pvalue)), size=0.3, shape= 21) +   # Show all points
theme_cowplot(12, font_size= 12) + #theme_minimal_hgrid(12, rel_small = -1) + 
scale_colour_gradientn(colors= colors_3) +
scale_fill_gradientn(colours= colors_3) +
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
    theme( legend.position="none") +
         xlab('Chromosome') +
    ylab('Z-score') +
	ylim(c(-5, 5)) +
geom_hline(yintercept= 0, size= 0.5, colour= 'black') #+
geom_hline(yintercept= -log10(HC), size= 0.5, linetype= 2, colour= '#878787') +
#annotate(geom="text", x= Inf, y= HC - 0.5, label= 'bold("High confidence")', color="black", vjust= 1, hjust= 1, parse= TRUE, size= 4)
geom_text_repel(data= genof, aes(x= BPcum, y= -log10(pvalue), label= gene), size= 3,  hjust = 0.5, force= 1, vjust= 1, colour= 'black')

save_plot(file= snakemake@output[[1]], plot= moms, base_width=297, base_height=210, units="mm", device= cairo_ps)

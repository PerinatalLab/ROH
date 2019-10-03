library(data.table)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(extrafont)

loadfonts()

d= fread(snakemake@input[[1]])
d$cf= 'HC'
d1= fread(snakemake@input[[2]])
d1$cf= 'LC'
d2= fread(snakemake@input[[3]])
d2$cf= 'NC'
d= rbind(d,d1,d2)
d= separate(d, segment, into= c('chr', 'cM1', 'cM2'), sep= ':')
d$cM1= as.numeric(d$cM1)
d$cM2= as.numeric(d$cM2)
d$chr= as.numeric(d$chr)

df= d

df$mcM= (df$cM1 + df$cM2) / 2

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

moms=  ggplot(don) +
    geom_point(aes(x=BPcum, y= z_meta, colour= as.factor(cf)), size=0.1) +   # Show all points
scale_colour_viridis_d(option='D', direction = 1) +
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
    #scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
#theme_bw() +
theme_bw() + # Custom the theme
    theme( text = element_text(size= 14),
           legend.position="none",
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()) +
	 xlab('Chromosome') +
    ylab('Z-score') +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))  +
geom_hline(yintercept= 0, size= 0.5, linetype= 2)


d= fread(snakemake@input[[4]])
d$cf= 'HC'
d1= fread(snakemake@input[[5]])
d1$cf= 'LC'
d2= fread(snakemake@input[[6]])
d2$cf= 'NC'
d= rbind(d,d1,d2)
d= separate(d, segment, into= c('chr', 'cM1', 'cM2'), sep= ':')
d$cM1= as.numeric(d$cM1)
d$cM2= as.numeric(d$cM2)
d$chr= as.numeric(d$chr)

df = d

df$mcM= (df$cM1 + df$cM2) / 2

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

fets=   ggplot(don) +
    geom_point(aes(x=BPcum, y= z_meta, colour= as.factor(cf)), size=0.1) +   # Show all points
scale_colour_viridis_d(option='D', direction = 1) +
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
    #scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
#theme_bw() +
theme_bw() + # Custom the theme
    theme( text = element_text(size= 14),
           legend.position="none",
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()) +
	 xlab('Chromosome') +
    ylab('Z-score') +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))  +
geom_hline(yintercept= 0, size= 0.5, linetype= 2)



#gg= arrangeGrob(moms,fets, ncol= 1)

ggsave(file= snakemake@output[[1]], plot= moms, dpi= 'retina', width= 18, height= 7, units= 'cm')
ggsave(file= snakemake@output[[2]], plot= fets, dpi= 'retina', width= 18, height= 7, units= 'cm')

library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(data.table)
library(tidyr)

colors_3= c('#FFBD01', '#00B25D', '#9C02A7')

x= fread(snakemake@input[[1]], h=T)

#x= separate(x, segment, into= c('cM1', 'cM2'), sep= ':')
#x= mutate(x, cM1= as.numeric(cM1), cM2= as.numeric(cM2))

x$mcM= (x$cM1 + x$cM2) / 2

x= x %>% filter(!is.na(chr))
  don1 <- x %>%
    group_by(chr)      %>%
    summarise(chr_len= max(mcM)) %>%
    mutate(tot= cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(x, ., by= 'chr') %>%
    arrange(chr, mcM) %>% # Add a cumulative position of each SNP
    mutate( BPcum=mcM+tot)

  axisdf = don1 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('chr', 'center')


mom= ggplot(don1) +
    geom_point(aes(x=BPcum, y= freq, colour= as.factor(chr)), size=0.1) +   # Show all points
theme_cowplot(12, font_size= 12) + 
scale_color_manual(values = rep(c(colors_3[1], '#00CCB5'), 23)) +
#scale_colour_viridis_d(option='D', direction = -1) + 
scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
    theme( text = element_text(size= 12),
           legend.position="none") +
         xlab('Chromosome') +
    ylab('Frequency') +
ggtitle('Mothers')

#return(mom)

x= fread(snakemake@input[[2]], h=T)

#x= separate(x, segment, into= c('cM1', 'cM2'), sep= ':')
#x= mutate(x, cM1= as.numeric(cM1), cM2= as.numeric(cM2))

x$mcM= (x$cM1 + x$cM2) / 2

x= x %>% filter(!is.na(chr))
  don1 <- x %>%
    group_by(chr)      %>%
    summarise(chr_len= max(mcM)) %>%
    mutate(tot= cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(x, ., by= 'chr') %>%
    arrange(chr, mcM) %>% # Add a cumulative position of each SNP
    mutate( BPcum=mcM+tot)

  axisdf = don1 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('chr', 'center')


dad= ggplot(don1) +
    geom_point(aes(x=BPcum, y= freq, colour= as.factor(chr)), size=0.1) +   # Show all points
theme_cowplot(12, font_size= 12) + 
scale_color_manual(values = rep(c('#EE3377','#DC125C'), 23)) +
#scale_colour_viridis_d(option='D', direction = -1) + 
scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
    theme( text = element_text(size= 12),
           legend.position="none") +
         xlab('Chromosome') +
    ylab('Frequency') +
ggtitle('Fathers')
#return(dad)


x= fread(snakemake@input[[3]], h=T)

#x= separate(x, segment, into= c('cM1', 'cM2'), sep= ':')
#x= mutate(x, cM1= as.numeric(cM1), cM2= as.numeric(cM2))

x$mcM= (x$cM1 + x$cM2) / 2

x= x %>% filter(!is.na(chr))
  don1 <- x %>%
    group_by(chr)      %>%
    summarise(chr_len= max(mcM)) %>%
    mutate(tot= cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(x, ., by= 'chr') %>%
    arrange(chr, mcM) %>% # Add a cumulative position of each SNP
    mutate( BPcum=mcM+tot)

  axisdf = don1 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('chr', 'center')

fet= ggplot(don1) +
    geom_point(aes(x=BPcum, y= freq, colour= as.factor(chr)), size=0.1) +   # Show all points
theme_cowplot(12, font_size= 12) +
scale_color_manual(values = rep(c('#c39a15','#E7B924'), 23)) +
#scale_colour_viridis_d(option='D', direction = -1) + 
scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
    theme( text = element_text(size= 12),
           legend.position="none") +
         xlab('Chromosome') +
    ylab('Frequency') +
ggtitle('Offspring')

plot1= plot_grid(mom, dad, fet, ncol= 1, labels= 'AUTO', align= 'v')

save_plot(snakemake@output[[1]], plot= plot1, base_width=297, base_height=210, units="mm", device= cairo_ps)

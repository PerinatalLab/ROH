library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)


cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

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

p1= ggplot(x, aes(x= as.numeric(as.character(SNP)), y= R2, colour= as.factor(pruning))) +
geom_point(size= 1) +
geom_line() +
scale_colour_viridis_d(name = "Pruning", labels = c("None", "0.9", "0.5", '0.1')) +
facet_grid(het~bp_cm, labeller = labeller(het= as_labeller(c('0'= 'No het. allowed', '1'= 'One het. allowed')), bp_cm= as_labeller(c('BP'='Physical distance', 'cM'='Genetic distance')))) +
 theme_bw(base_size=9, base_family = "Source Sans Pro") +
 theme(panel.border = element_blank(), 
 panel.grid.major = element_blank(), 
 panel.grid.minor = element_blank()) +
 theme(legend.position = "bottom") +
  theme(plot.title=element_text(size = 14)) +
  theme(axis.text.y=element_text(size=9)) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(axis.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9)) +
  theme(legend.text=element_text(size=9)) +
  theme(strip.text= element_text(size = 9)) +
	xlab('Number of genetic variants included') +
	ylab(expression(paste(R^2, ' between offspring ROH and parental genetic relatedness'))) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))+
ggtitle('A')

##### Figure 1.

bp= filter(d, bp_cm== 'BP')
cm= filter(d, bp_cm== 'cM')

colnames(bp)[1] <- "R2_bp"
colnames(cm)[1] <- "R2_cm"
df= inner_join(bp, cm, by= c('SNP', 'het', 'cohort', 'pruning'))

p2= ggplot(df, aes(x= R2_bp, y= R2_cm, colour= cohort)) +
geom_point(size= 1) +
geom_rug(col="grey",alpha=0.1, size=1.5)+
geom_abline(intercept = 0, slope = 1) +
scale_colour_viridis_d(name = "Sub-cohorts", labels = c("Cohort1", "Cohort2", "Cohort3", 'Cohort4', 'Cohort5', 'Cohort6')) +
 theme_bw(base_size=9, base_family = "Source Sans Pro") +
 theme(panel.border = element_blank(), 
 panel.grid.major = element_blank(), 
 panel.grid.minor = element_blank()) +
 theme(legend.position = "right") +
  theme(plot.title=element_text(size = 14)) +
  theme(axis.text.y=element_text(size=9)) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(axis.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9)) +
  theme(legend.text=element_text(size=9)) +
  theme(strip.text= element_text(size = 9)) +
	xlab(expression(paste(R^2, ' obtained using physical distance'))) +
	ylab(expression(paste(R^2, ' obtained using genetic distance'))) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2)) +
ggtitle('B')


lay <- rbind(c(1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1),
             c(NA,3,3,3,3,3,3,NA))


gg= arrangeGrob(p1, p2, ncol= 1)
ggsave(snakemake@output[[1]], gg, dpi= 'retina')

#grid.arrange(p1, p2, layout_matrix = lay)



#### Supp figure

s1= ggplot(filter(d, bp_cm== 'cM'), aes(x= as.numeric(as.character(SNP)), y= R2, colour= as.factor(pruning))) +
geom_point(size= 1) +
geom_line() +
scale_x_continuous(breaks=c(0, 200, 400)) +
scale_colour_viridis_d(name = "Pruning", labels = c("None", "0.9", "0.5", '0.1')) +
facet_grid(het~cohort, labeller = labeller(het= as_labeller(c('0'= 'No het. allowed', '1'= 'One het. allowed')), 
cohort= as_labeller(c('harvestm12'='Cohort1', 'harvestm24'='Cohort2', 'rotterdam1'='Cohort3', 'rotterdam2'='Cohort4', 'normentfeb'='Cohort5', 'normentmay'='Cohort6')))) +
 theme_bw(base_size=9, base_family = "Source Sans Pro") +
 theme(panel.border = element_blank(), 
 panel.grid.major = element_blank(), 
 panel.grid.minor = element_blank()) +
 theme(legend.position = "bottom") +
  theme(plot.title=element_text(size = 14)) +
  theme(axis.text.y=element_text(size=9)) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(axis.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9)) +
  theme(legend.text=element_text(size=9)) +
  theme(strip.text= element_text(size = 9)) +
	xlab('Number of genetic variants included') +
	ylab(expression(paste(R^2, ' between offspring ROH and parental genetic relatedness'))) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))


#ggsave(snakemake@output[[1]], grid.draw(g), device= 'eps', dpi= 'retina', width= 12, height= 8, units= 'cm')
ggsave(snakemake@output[[2]], plot= s1, device= 'eps', dpi= 'retina', width= 12, height= 10, units= 'cm')


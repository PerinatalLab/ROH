library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(circlize)

cohorts= c('harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay')

df_list= list()


input= unlist(snakemake@input)



for (coh in cohorts){
input_coh= unlist(input[grep(coh, input)])

trio= fread(unlist(input_coh[grep('trio', input_coh)]))


mom= fread(unlist(input_coh[grep('mfr_maternal', input_coh)]))
dad= fread(unlist(input_coh[grep('mfr_paternal', input_coh)]))
fet= fread(unlist(input_coh[grep('mfr_fetal', input_coh)]))

mom= select(mom, NSEG, FKB, KBAVG, IID)
dad= select(dad, NSEG, FKB, KBAVG, IID)
fet= select(fet, NSEG, FKB, KBAVG, IID)

names(mom)= paste0(names(mom), '_mom')
names(dad)= paste0(names(dad), '_dad')
names(fet)= paste0(names(fet), '_fet')

df= full_join(trio, mom, by= c('Mother'= 'IID_mom')) %>% full_join(., dad, by= c('Father'= 'IID_dad')) %>% full_join(., fet, by= c('Child'= 'IID_fet'))

df$cohort= coh

df= select(df, NSEG_mom, FKB_mom, KBAVG_mom, NSEG_dad, FKB_dad, KBAVG_dad, NSEG_fet, FKB_fet, KBAVG_fet, cohort)

print(summary(df))

df_list= c(df_list, list(df))

}

d= do.call('rbind', df_list)

matcor= melt(cor(d[,1:(ncol(d)-1)]))

matcor= filter(matcor, value> 0.1)

matcor$group1= ifelse(grepl('mom', matcor$from), 'mother', ifelse(grepl('dad', matcor$from), 'father', 'offspring'))
matcor$value1= ifelse(grepl('NSEG', matcor$from), 'NSEG', ifelse(grepl('FKB', matcor$from), 'Autozygosity', 'ROH average'))
matcor$group2= ifelse(grepl('mom', matcor$to), 'mother', ifelse(grepl('dad', matcor$to), 'father', 'offspring'))
matcor$value2= ifelse(grepl('NSEG', matcor$to), 'NSEG', ifelse(grepl('FKB', matcor$to), 'Autozygosity', 'ROH average'))
matcor$from= paste(matcor$value1, matcor$group1)
matcor$to= paste(matcor$value2, matcor$group2)
matcor= select(matcor, -group1, -group2, -value1, -value2)

matcor= filter(matcor, as.character(from)> as.character(to))

var_order= c("Autozygosity offspring", "ROH average offspring", "NSEG offspring", "Autozygosity mother", "ROH average mother", "NSEG mother", "Autozygosity father", "ROH average father", "NSEG father")

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- sort(rep(viridis(3, alpha = 1, begin = 0, end = 1, option = "D"), 3))

postscript(snakemake@output[[1]])
chordDiagram(
  x = matcor, 
order= var_order,
  grid.col = mycolor,
  transparency = 0.4,
  directional = 0,
#  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 2.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
      )

    # Add graduation on axis
   #" circos.axis(
    #  h = "top", 
    #  major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
    #  minor.ticks = 1, 
    #  major.tick.percentage = 0.5,
    #  labels.niceFacing = FALSE)
  }
)

dev.off()



x= d %>% gather(member, FKB, c('FKB_mom', 'FKB_dad', 'FKB_fet'))


p= ggplot(x, aes(sample = FKB, colour = factor(member))) +
  stat_qq(size= 1) +
  stat_qq_line() +
scale_colour_viridis_d(name = "Family member", labels = c("Father", "Offspring", "Mother")) +
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
	xlab('Theoretical quantile') +
	ylab('Estimated autozygosity') +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))

ggsave(snakemake@output[[2]])

library(data.table)
library(dplyr)
library(mclust)

d=fread(snakemake@input[[1]], h=T)

fit= Mclust(d$KB, 3)
d$ROH_class= fit$classification

x= group_by(d, ROH_class) %>% summarize(min_dist= min(KB, na.rm=T), max_dist= max(KB, na.rm=T))

x$min_distance= (shift(x$max_dist, type='lag') + x$min_dist)/ 2

x= filter(x, ROH_class>1) %>% select(ROH_class, min_distance)


write.table(filter(d, ROH_class> 2), snakemake@output[[1]], row.names= FALSE, col.names=T, sep= '\t', quote= F)
write.table(x, snakemake@output[[2]], row.names= FALSE, col.names=T, sep= '\t', quote=F)


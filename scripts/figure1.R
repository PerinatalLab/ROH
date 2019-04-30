library(data.table)
library(dplyr)
library(ggplot2)


pca= fread(snakemake@input[[1]])

fam_ibd= fread(snakemake@input[[2]])
names(fam_ibd)= c('FID', 'IID', 'x1','x2', 'x3','x4')

fam_roh= fread(snakemake@input[[3]])
names(fam_roh)= c('FID', 'IID', 'x1','x2', 'x3','x4')

ibd= fread(snakemake@input[[4]])
roh= fread(snakemake@input[[5]])

if (grepl('harvest', snakemake@input[[1]])) {
names(pca)= c('IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
} else {
names(pca)= c('FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
roh = rename(roh, SentrixID_1= SentrixID)
}




d= full_join(roh, ibd, by= c('SentrixID_1'= 'Child'))

d= filter(d, Mother %in% pca$IID,
		Father %in% pca$IID,
		Mother %in% fam_ibd$IID,
		SentrixID_1 %in% fam_roh$IID,
		Father %in% fam_ibd$IID,
		Mother %in% fam_ibd$IID)



d$cM= ifelse(is.na(d$cM), 0, d$cM)
d$FKB= ifelse(is.na(d$FKB), 0, d$FKB)


f1= ggplot(d, aes(FKB, cM)) +
geom_point() +
#scale_colour_gradient(low="#a1dab4", high="#41b6c4") +
theme_bw(base_family = "Source Sans Pro") +
theme( text = element_text(size= 14),
           legend.position="none",
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()
    ) + xlab('Parental relatedness, cM') +
    ylab('Offspring FROH, %') +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2))



ggsave(plot= f1, snakemake@output[[1]], device= 'eps', units= 'cm', dpi= 'print', width= 5, height= 5)



library('data.table')
library('ggplot2')
library('dplyr')
library('tidyr')
library('ggrepel')
library('cowplot')
library('survminer')
library(flexsurv)

colors_2= c('#9C02A7', '#00B25D')

d= fread(snakemake@input[[1]])

pvals= fread(snakemake@input[[2]])

if (nrow(pvals)==0) {
p1= ggplot() + theme_void()

save_plot(snakemake@output[[1]], plot= p1, base_width=297, base_height=210, units="mm", device= cairo_ps)

} else {


d= filter(d, paste(CHR, segment, sep=':') %in% pvals$segment)

rownames(d)= paste(d$CHR, d$segment, sep= ':')
cols= paste0('X', rownames(d))
cols= gsub(':', '_', cols)
df= as.data.frame(t(d[, 3:ncol(d)]))
ids= names(d[3:ncol(d)])
names(df)= cols

pvals= pvals[match(sub('X', '', cols), gsub(':', '_', pvals$segment)),]
pvals= separate(pvals, segment, into= c('chr', 'cM1', 'cM2'), sep= ':')
labels= paste(pvals$chr, pvals$pos1, pvals$pos2, sep=':')
df= as.data.frame(lapply(df, factor))
df$IID= ids

pheno= fread(snakemake@input[[3]])
df= inner_join(df, pheno, by= 'IID')
df= df[!duplicated(df$IID),]

surv_list= list()
for (i in 1:length(cols)) {
    
AFT_formula= as.formula(paste("Surv(SVLEN_UL_DG, spont) ~ ", paste0("`",cols[i],"`")))
fit_aft <- flexsurvreg(AFT_formula , dist = "weibull", data= df)
   
surv_covs <- summary(fit_aft, type = "survival", tidy = TRUE)
names(surv_covs)= c('time', 'est', 'lcl', 'ucl', 'ROH')
surv_covs$segment= labels[i]
surv_list[[i]] = surv_covs
}
  
x= do.call('rbind', surv_list)

x= filter(x, segment== cols[1])
 
p1= ggplot(data= x, aes(x = time, y = est, colour =  ROH )) + 
theme(legend.position = c(0, 225), 
legend.background=element_blank(),
legend.key=element_blank(),
legend.title=element_blank(),
strip.text = element_text(size= 12)) +  
geom_line() +
  facet_wrap(~segment) +
  scale_colour_manual(values= colors_2, breaks=c('0', '1'), labels=c('No autozygous segment', 'Autozygous segment')) +
  theme_cowplot(12, font_size= 12) + 
  xlab("Time to spontaneous delivery, days") + 
  ylab("Survival probability") +
  geom_ribbon(data= x, aes(ymin= lcl, ymax= ucl, fill= ROH), linetype=2, alpha=0.1, size= 0.4, show.legend= FALSE) +
  xlim(c(min(x$time), 308)) +
  geom_vline(xintercept= 260, linetype= 2, size= 0.4, colour= 'black', alpha= 0.6) +
theme(strip.background = element_blank()) 


save_plot(snakemake@output[[1]], plot= p1, base_width=297, base_height=210, units="mm", device= cairo_ps)
}

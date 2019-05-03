library(data.table)
library(dplyr)
library(tidyr)

options(stringsAsFactors=F)

fam= fread(snakemake@input[[1]])
names(fam)= c('FID', 'IID', 'x1', 'x2', 'x3', 'x4')

kob= fread(snakemake@input[[2]])
kob= gather(kob, ROLE, Retrieval, retrievalDetail_ID_B_315_926_930:retrievalDetail_ID_F_315_926_930)
kob$ROLE= ifelse(grepl('_B_', kob$ROLE), 'Child', ifelse(grepl('_M_', kob$ROLE), 'Mother', 'Father'))


rec= fread(snakemake@input[[3]])
names(rec)= c('false_FID', 'Retrieval', 'FID2', 'IID')

flag= fread(snakemake@input[[4]])
flag= filter(flag, phenoOK== TRUE, coreLMM== TRUE, genotypesOK== TRUE)

d= inner_join(fam, rec, by= 'IID')
d= inner_join(d, kob, by= 'Retrieval')

d= filter(d, FID== FID2)

df= select(d, PREG_ID_315, IID, x1, x2) 

names(df)= c('PREG_ID_315', 'Child', 'Father', 'Mother')
df= filter(df, !is.na(Child), !is.na(Father), !is.na(Mother))

write.table(df, snakemake@output[[1]], col.names=T, row.names=F, sep= '\t')

d= d %>% select(FID, IID, PREG_ID_315, ROLE)
names(d)= c('postFID', 'SentrixID', 'PREG_ID_315', 'Role')
write.table(d, snakemake@output[[2]], col.names= T, row.names=F, sep= '\t')

mfr= fread(snakemake@input[[5]])

write.table(mfr, snakemake@output[[3]], col.names=T, row.names=F, sep= '\t')

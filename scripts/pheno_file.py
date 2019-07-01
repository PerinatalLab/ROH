#!/usr/bin/python3

import pandas as pd
import numpy as np

wild= snakemake.wildcards.cohort

def pheno_harvest():
	d= pd.read_csv(snakemake.input[1], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[2], sep= ';', header= 0)
	link= pd.read_csv(snakemake.input[3], sep= ';', header= 0)
	pca= pd.read_csv(snakemake.input[4], sep= ' ', header= None, names= ['SentrixID_1', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'])
	bim= [line.strip() for line in open(snakemake.input[7], 'r')]
	bim= "".join(bim[1])
	bim= pd.read_csv(bim, sep= '\t', header= None, names= ['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	
	bp= bim.groupby(['chr'])['pos'].diff(1).sum() / 1000000
	
	
#	link= link.loc[link.PREG_ID_1724.str.isnumeric(), :]
#	link['PREG_ID_1724']= link['PREG_ID_1724'].astype(str).astype(int)
	link.dropna(subset= ['PREG_ID_1724'], inplace= True)
	mfr= pd.merge(mfr, link, on= ['PREG_ID_1724'], how= 'inner')
	mfr= pd.merge(mfr, pca, on= ['SentrixID_1'], how= 'inner')
	
	d= pd.merge(d, mfr, left_on= ['IID'], right_on= ['SentrixID_1'], how= 'inner')
	d['KB']= d['KB'] / 1000000 * 1000
	d[['FKB', 'FKBAVG', 'FNSEG']]= d[['KB', 'KBAVG', 'NSEG']].divide(bp, axis=1)

	fam= pd.read_csv(snakemake.input[6], sep=' ', header= None)
	fam.columns= ['FID','IID','m','f','sex','pheno']
	
	d= d[(d['FLERFODSEL']==0)]
	d= d[(d['DODKAT']<6) | (d['DODKAT']>10)]
	d= d[(d['SVLEN_UL_DG']< '308')]
	d= d[(d['SVLEN_UL_DG']> '154')]
	d['SVLEN_UL_DG']= d['SVLEN_UL_DG'].replace(' ', np.nan)
	d['SVLEN_UL_DG']= d.SVLEN_UL_DG.astype(str).astype(int)
	d= d.dropna(subset= ['SVLEN_UL_DG'])
	d['IVF']= d['IVF'].replace(' ', np.nan)
	d= d[(pd.isnull(d['IVF']))] 
	d= d[(d.ABRUPTIOP==0)]
	d= d[(d.PLACENTA_PREVIA==0) ]
	d= d[(d.FOSTERV_POLYHYDRAMNION==0)]
	d= d[(d.C00_MALF_ALL==0)]
	d['BATCH']= np.where(d.SentrixID_1.isin(fam.IID), 0, 1)
	d= d.sample(frac=1)
	flag= pd.read_csv(snakemake.input[8], sep= '\t', header= 0)
	flag= flag[(flag['genotypesOK']== True) & (flag['phenotypesOK']== True) & (flag['coreOK']== True)]
	d= d.loc[d.IID.isin(flag.IID), :]
	d.drop_duplicates(subset= ['PREG_ID_1724'], keep= 'first', inplace= True)
	return d

def pheno_rotterdam():
	d= pd.read_csv(snakemake.input[1], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
	link= pd.read_csv(snakemake.input[3], delim_whitespace= True, header= 0)
	pca= pd.read_csv(snakemake.input[4], delim_whitespace= True, header= 0)
	pca.columns= ['FID', 'SentrixID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
	bim= [line.strip() for line in open(snakemake.input[7], 'r')]
	bim= "".join(bim[1])
	bim= pd.read_csv(bim, sep= '\t', header= None, names=['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	bp= bim.groupby(['chr'])['pos'].diff(1).sum() / 1000000
	mfr= pd.merge(mfr, link, on= ['PREG_ID_315'], how= 'inner')
	mfr= pd.merge(mfr, pca, on= ['SentrixID'], how= 'inner')
	
	d= pd.merge(d, mfr, left_on= ['IID'], right_on= ['SentrixID'], how= 'inner')
	d['KB']= d['KB'] / 1000000 * 1000
	d[['FKB', 'FKBAVG', 'FNSEG']]= d[['KB', 'KBAVG', 'NSEG']].divide(bp, axis=1)
	d= d[(d['FLERFODSEL']=='Enkeltfødsel')]
	d= d[d['DODKAT'].str.contains('Levendefødt')]
	d= d[(d['SVLEN_UL_DG']< 308)]
	d= d[(d['SVLEN_UL_DG']> 154)]
	d= d.dropna(subset= ['SVLEN_UL_DG'])
	d['IVF']= d['IVF'].replace(' ', np.nan)
	d= d[(pd.isnull(d['IVF']))]
	d= d[(d.ABRUPTIOP=='Nei')]
	d= d[(d.PLACENTA_PREVIA=='Nei') ]
	d= d[(d.FOSTERV_POLYHYDRAMNION=='Nei')]
	d= d[(d.C00_MALF_ALL=='Nei')]
	flag= pd.read_csv(snakemake.input[8], sep= '\t', header= 0)
	flag= flag[(flag['genotypesOK']== True) & (flag['phenoOK']== True) & (flag['coreLMM']== True)]
	d= d.loc[d.IID.isin(flag.IID), :]
	d.drop_duplicates(subset= ['PREG_ID_315'], keep= 'first', inplace= True)
	return d

### Relatedness

def selectUnrelated(df, x):
	kin= pd.read_csv(snakemake.input[5], sep= '\t')
	kin= kin.loc[kin.KINSHIP > 0.0884, :]
	kin= kin.loc[kin.ID1.isin(x.values)]
	kin= kin.loc[kin.ID2.isin(x.values)]
	kin= kin.loc[:, ['ID1','ID2','KINSHIP']]
	kin_temp= kin.copy()
	kin_temp.columns= ['ID2', 'ID1', 'KINSHIP']
	kin_temp= kin_temp.append(kin)
	kin_temp['n']= kin_temp.groupby('ID1')['ID1'].transform('count')
	kin_temp['nn']= kin_temp.groupby('ID2')['ID2'].transform('count')
	kin_temp.sort_values(by=['n', 'nn'], inplace= True)
	to_keep= list()
	for i in range(0, len(kin_temp.index)):
		if kin_temp.iloc[i, 0] in kin_temp.iloc[0:i, 1].values:
			kin_temp.iloc[i, 1]= "X"
		else:
			to_keep.append(kin_temp.iloc[i, 0])
	to_remove= [i for i in kin_temp.ID1 if i not in to_keep]
	to_remove= list(set(to_remove))
	remove= pd.DataFrame({'FID': to_remove})
	remove['IID']= remove.FID
	return remove

if 'harvest' in wild:
	d= pheno_harvest()
	remove= selectUnrelated(d, d.SentrixID_1)
	d= d[~d.SentrixID_1.isin(remove)]
if 'harvest' not in wild:
	d= pheno_rotterdam()
	remove= selectUnrelated(d, d.SentrixID)
	d= d[~d.SentrixID.isin(remove)]

d.to_csv(snakemake.output[0], sep= '\t', index= False)



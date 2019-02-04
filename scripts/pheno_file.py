#!/usr/bin/python3

import pandas as pd
import numpy as np

wild= snakemake.wildcards.cohort

def pheno_harvest():
	d12= pd.read_csv(snakemake.input[1], delim_whitespace= True)
	d24= pd.read_csv(snakemake.input[3], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[4], sep= ';', header= 0)
	link= pd.read_csv(snakemake.input[5], sep= ';', header= 0)
	pca= pd.read_csv(snakemake.input[6], sep= ' ', header= None, names= ['SentrixID_1', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'])
	bim12= pd.read_csv(snakemake.input[7], sep= '\t', header= None, names=['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	bim24= pd.read_csv(snakemake.input[8], sep= '\t', header= None, names= ['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	
	bp12= bim12.groupby(['chr'])['pos'].diff(1).sum()
	bp24= bim24.groupby(['chr'])['pos'].diff(1).sum()
	
	pca.columns= ['SentrixID_1', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
	
	link= link.loc[link.PREG_ID_1724.str.isnumeric(), :]
	link['PREG_ID_1724']= link['PREG_ID_1724'].astype(str).astype(int)
	
	mfr= pd.merge(mfr, link, on= ['PREG_ID_1724'], how= 'inner')
	mfr= pd.merge(mfr, pca, on= ['SentrixID_1'], how= 'inner')
	
	d12= pd.merge(d12, mfr, left_on= ['IID'], right_on= ['SentrixID_1'], how= 'inner')
	d24= pd.merge(d24, mfr, left_on= ['IID'], right_on= ['SentrixID_1'], how= 'inner')
	
	d12[['FKB', 'FKBAVG', 'FNSEG']]= d12[['KB', 'KBAVG', 'NSEG']].divide(bp12, axis=1)
	d24[['FKB', 'FKBAVG', 'FNSEG']]= d24[['KB', 'KBAVG', 'NSEG']].divide(bp24, axis=1)
	
	d= pd.concat([d12, d24])
	
	fam12= pd.read_csv(snakemake.input[14], sep='\t', header= None)
	
	runs12= pd.read_csv(snakemake.input[0], delim_whitespace= True)
	runs24= pd.read_csv(snakemake.input[2], delim_whitespace= True)
	runs= pd.concat([runs12, runs24], axis= 0)
	runs['long']= np.where((runs.POS2 - runs.POS1) > 8500000, 1, 0)
	runs['bp']= runs.POS2 - runs.POS1
	x= runs.groupby(['IID']).sum()[['long']]
	x_l= runs.loc[runs.long==1, :].groupby(['IID']).sum()[['bp']]
	x_s= runs.loc[runs.long==0, :].groupby(['IID']).sum()[['bp']]
	x.columns= ['long_ROH_count']
	x_l.columns= ['KB_long']
	x_s.columns= ['KB_short']
	df= pd.concat([x, x_l, x_s], axis= 1)
	df['IID']= df.index
	df.fillna(0, inplace=True)
	d= pd.merge(d, df, left_on= ['SentrixID_1'], right_on= ['IID'])
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
	d['BATCH']= np.where(d.SentrixID_1.isin(fam12.IID), 0, 1)
	d= d.sample(frac=1)
	
	d.drop_duplicates(subset= ['PREG_ID_1724'], keep= 'first', inplace= True)
	return d

def pheno_rotterdam1():
	d= pd.read_csv(snakemake.input[11], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[12], sep= '\t', header= 0)
	link= pd.read_csv(snakemake.input[5], sep= ' ', header= 0)
	pca= pd.read_csv(snakemake.input[6], sep= ' ', header= None, names= ['SentrixID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'])
	bim= pd.read_csv(snakemake.input[13], sep= '\t', header= None, names=['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	bp= bim.groupby(['chr'])['pos'].diff(1).sum()
	mfr= pd.merge(mfr, link, on= ['PREG_ID_315'], how= 'inner')
	mfr= pd.merge(mfr, pca, on= ['SentrixID'], how= 'inner')
	
	d= pd.merge(d, mfr, left_on= ['IID'], right_on= ['SentrixID'], how= 'inner')
	d[['FKB', 'FKBAVG', 'FNSEG']]= d[['KB', 'KBAVG', 'NSEG']].divide(bp, axis=1)
	runs= pd.read_csv(snakemake.input[10], delim_whitespace= True)
	runs['long']= np.where((runs.POS2 - runs.POS1) > 8500000, 1, 0)
	runs['bp']= runs.POS2 - runs.POS1
	x= runs.groupby(['IID']).sum()[['long']]
	x_l= runs.loc[runs.long==1, :].groupby(['IID']).sum()[['bp']]
	x_s= runs.loc[runs.long==0, :].groupby(['IID']).sum()[['bp']]
	x.columns= ['long_ROH_count']
	x_l.columns= ['KB_long']
	x_s.columns= ['KB_short']
	
	df= pd.concat([x, x_l, x_s], axis= 1)
	df['IID']= df.index
	df.fillna(0, inplace=True)
	d= pd.merge(d, df, left_on= ['SentrixID'], right_on= ['IID'])
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
	d.drop_duplicates(subset= ['PREG_ID_315'], keep= 'first', inplace= True)
	return d


### Relatedness

def selectUnrelated(df, x):
	kin= pd.read_csv(snakemake.input[9], sep= '\t')
	kin= kin.loc[kin.KINSHIP > 0.0884, :]
	kin= kin.loc[kin.ID1.isin(x.values)]
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
	remove= pd.DataFrame(to_remove)
	remove.columns= ['FID']
	remove['IID']= remove.FID
	return remove

if wild == 'harvest':
	d= pheno_harvest()
	remove= selectUnrelated(d, d.SentrixID_1)
	d= d[~d.SentrixID_1.isin(remove)]
elif wild == 'rotterdam1':
	d= pheno_rotterdam1()
	remove= selectUnrelated(d, d.SentrixID)
	d= d[~d.SentrixID.isin(remove)]

d.to_csv(snakemake.output[0], sep= '\t', index= False)



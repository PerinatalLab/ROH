#!/usr/bin/python3

import pandas as pd
import numpy as np
from sklearn.preprocessing import scale

coh= snakemake.wildcards.cohort

def pheno_harvest():
	d= pd.read_csv(snakemake.input[1], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
	link= pd.read_csv(snakemake.input[3], sep= '\t', header= 0)
	bim= [line.strip() for line in open(snakemake.input[7], 'r')]
	bim= "".join(bim[1])
	bim= pd.read_csv(bim, sep= '\t', header= None, names= ['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	
	bp= bim.groupby(['chr'])['pos'].diff(1).sum() / 1000000
	
	
	link.dropna(subset= ['PREG_ID_1724'], inplace= True)
	mfr= pd.merge(mfr, link, on= ['PREG_ID_1724'], how= 'inner')
	d= pd.merge(d, mfr, left_on= ['IID'], right_on= ['SentrixID_1'], how= 'inner')
	d['KB']= d['KB'] / 1000000 * 1000
	d[['FKB', 'FKBAVG', 'FNSEG']]= d[['KB', 'KBAVG', 'NSEG']].divide(bp, axis=1)

	fam= pd.read_csv(snakemake.input[6], sep=' ', header= None)
	fam.columns= ['FID','IID','m','f','sex','pheno']
	
	d= d[(d['FLERFODSEL']==0)]
	d= d[(d['DODKAT']<6) | (d['DODKAT']>10)]
	d= d[(d.ABRUPTIOP==0)]
	d= d[(d.PLACENTA_PREVIA==0) ]
	d['SVLEN_UL_DG']= np.where(d.SVLEN_UL_DG.isnull(), d.SVLEN_SM_DG, d.SVLEN_UL_DG)
	d= d[(d['SVLEN_UL_DG']< 308)]
	d= d[(d['SVLEN_UL_DG']> 154)]
	d= d.dropna(subset= ['SVLEN_UL_DG'])
	d= d[(pd.isnull(d['IVF']))]
	d= d[(d.FOSTERV_POLYHYDRAMNION==0)]
	d= d[(d.FOSTERV_OLIGOHYDRAMNION== 0)]
	d= d[(d.C00_MALF_ALL==0)]
	d= d[d.DIABETES_MELLITUS.isna()]
	d= d[(d.HYPERTENSJON_KRONISK==0)]
	d= d[(d.HYPERTENSJON_ALENE==0)]
	d= d[d.PREEKL.isna()]
	d= d.sample(frac=1)
	flag= pd.read_csv(snakemake.input[8], sep= '\t', header= 0)
	flag= flag[(flag['genotypesOK']== True) & (flag['phenotypesOK']== True)]
	d= d.loc[d.IID.isin(flag.IID), :]
	d.drop_duplicates(subset= ['PREG_ID_1724'], keep= 'first', inplace= True)
	pca_out= [line.strip() for line in open(snakemake.input[9], 'rt')]
	d['spont']= np.where((d.FSTART==1) | (((d.KSNITT.isnull()) | (d.KSNITT>1)) & ((d.KSNITT_PLANLAGT.isnull()) | (d.KSNITT_PLANLAGT==1)) & (d.INDUKSJON_PROSTAGLANDIN==0) & (d.INDUKSJON_ANNET==0) & (d.INDUKSJON_OXYTOCIN==0) & (d.INDUKSJON_AMNIOTOMI==0)) | (~d.VANNAVGANG.isnull()) , 1, 0)
	d['PARITY0']= np.where(d.PARITET_5==0, 1, 0)
	sUPD= [line.strip() for line in open(snakemake.input[10], 'rt')]
	d= d.loc[~d.IID.isin(pca_out), :]
	d= d.loc[~d.IID.isin(sUPD), :]
	d['SVLEN_UL_DG']= np.where(d.FKB< 0.08, d.SVLEN_UL_DG, np.nan)
	d= d.loc[d.VEKT> 1500, :]
	d['cohort']= coh
	d['PREG_ID']= d.PREG_ID_1724
	if snakemake.wildcards.sample== 'maternal':
		d['ALDER']= d.MORS_ALDER
	if snakemake.wildcards.sample== 'paternal':
		d['ALDER']= d.FARS_ALDER_KAT_K8
	if snakemake.wildcards.sample== 'fetal':
		d['ALDER']= d.FAAR
	d.drop_duplicates(subset= ['IID'], keep= 'first', inplace= True)
	d= d[['IID', 'PREG_ID', 'SVLEN_UL_DG', 'spont', 'ALDER', 'cohort', 'FKB', 'NSEG', 'KBAVG', 'PARITY0']]
	return d

def pheno_rotterdam():
	d= pd.read_csv(snakemake.input[1], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
	link= pd.read_csv(snakemake.input[3], delim_whitespace= True, header= 0)
	bim= [line.strip() for line in open(snakemake.input[7], 'r')]
	bim= "".join(bim[1])
	bim= pd.read_csv(bim, sep= '\t', header= None, names=['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	bp= bim.groupby(['chr'])['pos'].diff(1).sum() / 1000000
	mfr= pd.merge(mfr, link, on= ['PREG_ID_315'], how= 'inner')
	
	d= pd.merge(d, mfr, left_on= ['IID'], right_on= ['SentrixID'], how= 'inner')
	d['KB']= d['KB'] / 1000000 * 1000
	d[['FKB', 'FKBAVG', 'FNSEG']]= d[['KB', 'KBAVG', 'NSEG']].divide(bp, axis=1)
	d= d[(d['FLERFODSEL']=='Enkeltfødsel')]
	d= d[d['DODKAT'].str.contains('Levendefødt')]
	d= d[(d.ABRUPTIOP=='Nei')]
	d= d[(d.PLACENTA_PREVIA=='Nei') ]
	d['SVLEN_UL_DG']= np.where(d.SVLEN_UL_DG.isnull(), d.SVLEN_SM_DG, d.SVLEN_UL_DG)
	d= d[(d['SVLEN_UL_DG']< 308)]
	d= d[(d['SVLEN_UL_DG']> 154)]
	d= d.dropna(subset= ['SVLEN_UL_DG'])
	d['IVF']= d['IVF'].replace(' ', np.nan)
	d= d[(pd.isnull(d['IVF']))]
	d= d[(d.FOSTERV_POLYHYDRAMNION=='Nei')]
	d= d[(d.FOSTERV_OLIGOHYDRAMNION== 'Nei')]
	d= d[(d.C00_MALF_ALL=='Nei')]
	d= d[d.DIABETES_MELLITUS.isna()]
	d= d[(d.HYPERTENSJON_KRONISK=='Nei')]
	d= d[(d.HYPERTENSJON_ALENE=='Nei')]
	d= d[d.PREEKL.isna()]
	flag= pd.read_csv(snakemake.input[8], sep= '\t', header= 0)
	flag= flag[(flag['genotypesOK']== True) & (flag['phenoOK']== True)]
	d= d.loc[d.IID.isin(flag.IID), :]
	d.drop_duplicates(subset= ['PREG_ID_315'], keep= 'first', inplace= True)
	pca_out= [line.strip() for line in open(snakemake.input[9], 'rt')]
	d= d.loc[~d.IID.isin(pca_out), :]
	d['spont']= np.where(((d.FSTART=='Spontan') | (d.FSTART== '')) | (((d.KSNITT=='') | (d.KSNITT== 'Uspesifisert') | (d.KSNITT== 'Akutt keisersnitt')) & (d.INDUKSJON_PROSTAGLANDIN=='Nei') & (d.INDUKSJON_ANNET=='Nei') & (d.INDUKSJON_OXYTOCIN=='Nei') & (d.INDUKSJON_AMNIOTOMI=='Nei')) | (~d.VANNAVGANG.isnull()), 1, 0)
	d['PARITY0']= np.where(d.PARITET_5=='0 (førstegangsfødende)', 1, 0)
	sUPD= [line.strip() for line in open(snakemake.input[10], 'rt')]
	d= d.loc[~d.IID.isin(sUPD), :]
	d['cohort']= coh
	d['PREG_ID']= d.PREG_ID_315
	d['MORS_ALDER']= d.FAAR - d.MOR_FAAR
	if snakemake.wildcards.sample== 'maternal':
		d['ALDER']= d.FAAR - d.MOR_FAAR
	if snakemake.wildcards.sample== 'paternal':
		d['ALDER']= pd.Categorical(d.FARS_ALDER_KAT_K8).codes + 1
		d['ALDER']= np.where(d.ALDER== 0, np.nan, d.ALDER)
	if snakemake.wildcards.sample== 'fetal':
		d['ALDER']= d.FAAR
	d['SVLEN_UL_DG']= np.where(d.FKB< 0.08, d.SVLEN_UL_DG, np.nan)
	d= d.loc[d.VEKT> 1500, :]
	d.drop_duplicates(subset= ['IID'], keep= 'first', inplace= True)
	d= d[['IID', 'PREG_ID', 'SVLEN_UL_DG', 'spont', 'ALDER', 'cohort', 'FKB', 'NSEG', 'KBAVG', 'PARITY0']]
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

if 'harvest' in coh:
	d= pheno_harvest()
#	remove= selectUnrelated(d, d.SentrixID_1)
#	d= d[~d.SentrixID_1.isin(remove)]
if 'harvest' not in coh:
	d= pheno_rotterdam()
#	remove= selectUnrelated(d, d.SentrixID)
#	d= d[~d.SentrixID.isin(remove)]

d.to_csv(snakemake.output[0], sep= '\t', index= False, header= True)



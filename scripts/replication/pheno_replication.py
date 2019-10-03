#!/usr/bin/python3

import pandas as pd
import numpy as np

wild= snakemake.wildcards.rep

def pheno():
	mfr= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
	link= pd.read_csv(snakemake.input[1], delim_whitespace= True, header= 0)
	pca= pd.read_csv(snakemake.input[2], delim_whitespace= True, header= 0)
	pca.columns= ['FID', 'SentrixID', 'NMISS_ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
	mfr= pd.merge(mfr, link, on= ['PREG_ID_315'], how= 'inner')
	d= pd.merge(mfr, pca, on= ['SentrixID'], how= 'inner')
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
	flag= pd.read_csv(snakemake.input[4], sep= '\t', header= 0)
	flag= flag[(flag['genotypesOK']== True) & (flag['phenoOK']== True) & (flag['coreLMM']== True)]
	d= d.loc[d.SentrixID.isin(flag.IID), :]
	d.drop_duplicates(subset= ['PREG_ID_315'], keep= 'first', inplace= True)
	pca_out= [line.strip() for line in open(snakemake.input[5], 'rt')]
	d= d.loc[~d.SentrixID.isin(pca_out), :]
	return d

### Relatedness

def selectUnrelated(df, x):
	kin= pd.read_csv(snakemake.input[3], sep= '\t')
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

d= pheno()
remove= selectUnrelated(d, d.SentrixID)
d= d[~d.SentrixID.isin(remove)]

d.to_csv(snakemake.output[0], sep= '\t', index= False)



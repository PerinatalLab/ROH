#!/usr/bin/python3

import pandas as pd
import numpy as np

coh= snakemake.wildcards.rep

def pheno():
	d= pd.read_csv(snakemake.input[1], delim_whitespace= True)
	mfr= pd.read_csv(snakemake.input[2], sep='\t', header= 0)
	link= pd.read_csv(snakemake.input[3], delim_whitespace= True, header= 0)
	bim= [line.strip() for line in open(snakemake.input[4], 'r')]
	bim= "".join(bim[1])
	bim= pd.read_csv(bim, sep= '\t', header= None, names=['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])
	bp= bim.groupby(['chr'])['pos'].diff(1).sum() / 1000000
	mfr= pd.merge(mfr, link, on= 'PREG_ID_315', how= 'inner')
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
	d= d[(d.FOSTERV_POLYHYDRAMNION=='Nei')]
	d= d[(d.FOSTERV_OLIGOHYDRAMNION== 'Nei')]
	d= d[(d.C00_MALF_ALL=='Nei')]
	d= d[d.DIABETES_MELLITUS.isna()]
        d= d[(d.HYPERTENSJON_KRONISK=='Nei')]
        d= d[(d.HYPERTENSJON_ALENE=='Nei')]
        d= d[d.PREEKL.isna()]
	flag= pd.read_csv(snakemake.input[5], sep= '\t', header= 0)
	flag= flag[(flag['genotypesOK']== True) & (flag['phenoOK']== True)]
	d= d.loc[d.SentrixID.isin(flag.IID), :]
	d.drop_duplicates(subset= ['PREG_ID_315'], keep= 'first', inplace= True)
	pca_out= [line.strip() for line in open(snakemake.input[6], 'rt')]
	d['spont']= np.where(((d.FSTART=='Spontan') | (d.FSTART== '')) | (((d.KSNITT=='') | (d.KSNITT== 'Uspesifisert') | (d.KSNITT== 'Akutt keisersnitt')) & (d.INDUKSJON_PROSTAGLANDIN=='Nei') & (d.INDUKSJON_ANNET=='Nei') & (d.INDUKSJON_OXYTOCIN=='Nei') & (d.INDUKSJON_AMNIOTOMI=='Nei')) | (~d.VANNAVGANG.isnull()), 1, 0)
	d['PARITY0']= np.where(d.PARITET_5=='0 (førstegangsfødende)', 1, 0)
	sUPD= [line.strip() for line in open(snakemake.input[7], 'rt')]
	d= d.loc[~d.IID.isin(sUPD), :]
	d['PREG_ID']= d.PREG_ID_315
	d['MORS_ALDER']= d.FAAR - d.MOR_FAAR
	if snakemake.wildcards.sample_rep== 'maternal':
		d['ALDER']= d.FAAR - d.MOR_FAAR
	if snakemake.wildcards.sample_rep== 'paternal':
		d['ALDER']= pd.Categorical(d.FARS_ALDER_KAT_K8).codes + 1
		d['ALDER']= np.where(d.ALDER== 0, np.nan, d.ALDER)
	if snakemake.wildcards.sample_rep== 'fetal':
		d['ALDER']= d.FAAR
	d['SVLEN_UL_DG']= np.where(d.FKB< 0.08, d.SVLEN_UL_DG, np.nan)
	d= d.loc[d.VEKT> 1500, :]
	d['cohort']= coh
	d= d.loc[~d.SentrixID.isin(pca_out), :]
	d.drop_duplicates(subset= ['IID'], keep= 'first', inplace= True)
	d= d[['IID', 'PREG_ID', 'SVLEN_UL_DG', 'spont', 'ALDER', 'cohort', 'FKB', 'NSEG', 'KBAVG', 'PARITY0']]
	return d

d= pheno()

d.to_csv(snakemake.output[0], sep= '\t', index= False)



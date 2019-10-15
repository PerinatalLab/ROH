import pandas as pd
import numpy as np

df_list= list()

if 'frequency' not in snakemake.input[0]:
	d= pd.read_csv(snakemake.input[0], sep= '\t', header= None, names= ['segment', 'n', 'beta', 'sd', 'pvalue', 'R', 'R_pvalue'])
	d[['chr', 'cM1', 'cM2']]= d['segment'].str.split(':', expand= True)
	d[['chr', 'cM1', 'cM2']]= d[['chr', 'cM1', 'cM2']].apply(lambda x: x.astype('float'))
	for infile in snakemake.input[1:]:
		df= pd.read_csv(infile, sep= '\t', header= None, names= ['segment', 'n', 'beta', 'sd', 'pvalue', 'R', 'R_pvalue'])
		df[['chr', 'cM1', 'cM2']]= df['segment'].str.split(':', expand= True)
		df[['chr', 'cM1', 'cM2']]= df[['chr', 'cM1', 'cM2']].apply(lambda x: x.astype('float'))
		df= df[['chr', 'cM1', 'cM2']]
		df_list.append(df)
	
	df= pd.concat(df_list)
	
	df_list= list()
	for CHR in set(d.chr):
		a= df.loc[df.chr== CHR, :]
		a= pd.concat([a.cM1, a.cM2])
		a= np.unique(a)
		a= np.sort(a)
		temp_d= d.loc[d.chr== CHR, :]
		for index, row in temp_d.iterrows():
			bh= row.cM2
			bl= row.cM1
			i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
			x= pd.DataFrame(a[i], columns= ['cM1']).dropna()
			x['cM2']= x.cM1.shift(-1)
			x.dropna(inplace= True)
			x['chr'], x['n'], x['beta'], x['sd'], x['pvalue'], x['R'], x['R_pvalue']= row.chr, row.n, row.beta, row.sd, row.pvalue, row.R, row.R_pvalue
			df_list.append(x.copy())

if 'frequency' in snakemake.input[0]:
	d= pd.read_csv(snakemake.input[0], sep= '\t', header= None, names= ['chr', 'segment', 'freq'])
	d[['cM1', 'cM2']]= d['segment'].str.split(':',expand=True)
	d[['cM1', 'cM2']]= d[['cM1', 'cM2']].apply(lambda x: x.astype('float'))
	df_list= list()
	for infile in snakemake.input[1:]:
                df= pd.read_csv(infile, sep= '\t', header= None, names= ['chr', 'segment', 'freq'])
                df[['cM1', 'cM2']]= df['segment'].str.split(':', expand= True)
                df[['cM1', 'cM2']]= df[['cM1', 'cM2']].apply(lambda x: x.astype('float'))
                df= df[['chr', 'cM1', 'cM2']]
                df_list.append(df)
	df= pd.concat(df_list)
	
	df_list= list()
	for CHR in set(d.chr):
		a= df.loc[df.chr== CHR, :]
		a= pd.concat([a.cM1, a.cM2])
		a= np.unique(a)
		a= np.sort(a)
		temp_d= d.loc[d.chr== CHR, :]
		for index, row in temp_d.iterrows():
			bh= row.cM2
			bl= row.cM1 
			i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
			x= pd.DataFrame(a[i], columns= ['cM1']).dropna()
			x['cM2']= x.cM1.shift(-1)
			x.dropna(inplace= True)
			x['chr'], x['freq']= row.chr, row.freq
			df_list.append(x.copy())

df= pd.concat(df_list)

df.to_csv(snakemake.output[0], header= True, sep= '\t', index= False)

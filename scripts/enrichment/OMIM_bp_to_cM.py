import pandas as pd
import numpy as np

d= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
d['chr']= pd.to_numeric(d.chr.str.replace('chr', ''), errors= 'coerce')
d.dropna(subset= ['chr'], inplace= True)

g= pd.read_csv(snakemake.input[1], sep= ' ', header= 0)
g= g[['chr', 'Genetic_Map(cM)', 'position']]
g.columns= ['chr', 'cM', 'pos']

df= pd.merge(d, g, left_on= ['chr', 'start'], right_on= ['chr', 'pos'], how= 'left')

df_miss= df.loc[df['cM'].isna(), :]
newdf= pd.DataFrame()

for CHR in set(df_miss.chr):
	df_temp= df_miss.loc[df_miss.chr== CHR, :]
	g_temp= g.loc[g.chr== CHR, :]
	df_temp['newX']= np.interp(df_temp['start'], g_temp['pos'], g_temp['cM'])
	newdf= newdf.append(df_temp)

newdf= newdf[['chr', 'start', 'newX']]
df= pd.merge(df, newdf, on= ['chr', 'start'], how= 'left')
df['cM1_OMIM']= np.where(df['cM'].isna(), df['newX'], df['cM'])
df['cM1_OMIM']= (df.cM1_OMIM*10**4).round() * 100
df= df[['chr', 'Phenotype', 'end', 'cM1_OMIM']]
df= pd.merge(df, g, left_on= ['chr', 'end'], right_on= ['chr', 'pos'], how= 'left')
df_miss= df.loc[df['cM'].isna(), :]
newdf= pd.DataFrame()
for CHR in set(df_miss.chr):
	df_temp= df_miss.loc[df_miss.chr== CHR, :]
	g_temp= g.loc[g.chr== CHR, :]
	df_temp['newX']= np.interp(df_temp['end'], g_temp['pos'], g_temp['cM'])
	newdf= newdf.append(df_temp)

newdf= newdf[['chr', 'end', 'newX']]
df= pd.merge(df, newdf, on= ['chr', 'end'], how= 'left')
df['cM2_OMIM']= np.where(df['cM'].isna(), df['newX'], df['cM'])
df['cM2_OMIM']= (df.cM2_OMIM*10**4).round() * 100
df= df[['chr', 'Phenotype', 'cM1_OMIM', 'cM2_OMIM']]
df.dropna(subset= ['Phenotype'], inplace= True)
df.drop_duplicates(keep= 'first', inplace= True)
rec= df[df.Phenotype.str.contains('recessive', na= False)]
dom= df[df.Phenotype.str.contains('dominant', na= False)]
rec.to_csv(snakemake.output[0], sep= '\t', header= True, index= False)
dom.to_csv(snakemake.output[1], sep= '\t', header= True, index= False)

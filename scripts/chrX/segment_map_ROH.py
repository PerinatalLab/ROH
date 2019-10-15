import pandas as pd
import numpy as np

d= pd.read_csv(snakemake.input[0], delim_whitespace= True, header= 0)

x= [line.strip() for line in open(snakemake.input[1], 'r')]
fam= [i for i in x if 'fam' in i]

fam= pd.read_csv("".join(fam), delim_whitespace=True, header= None)

d= d.loc[:, ['IID','CHR','POS1','POS2']]

fam.columns= ['FID', 'IID', 'X1', 'X2', 'X3', 'X4']

#CHR= snakemake.wildcards.CHR

#chr_df= d.loc[d.CHR== int(float(CHR)),:]

a= pd.concat([d.POS1, d.POS2])
a= np.unique(a)
a= np.sort(a)


df_list= list()

for id in set(d.IID):
	temp_df= d.loc[d.IID== id, :]	
	for index, row in temp_df.iterrows():
		bh= row.POS2
		bl= row.POS1
		i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
		x= pd.DataFrame(a[i], columns= ['start']).dropna()
		x['end']= x.start.shift(-1)
		x.dropna(inplace= True)
		x['IID']= id
		df_list.append(x.copy())

df= pd.concat(df_list)
df['Value']= 1
df['segment']= df['start'].astype(str) + ':' + df['end'].astype(str)
df= df.pivot(index= 'segment', columns= 'IID', values= 'Value')
df['segment']= df.index
df['CHR']= 23

cols = ['CHR', 'segment']  + [col for col in df if col not in 'CHR,segment']
df= df[cols]

o_ids= fam[~fam.IID.isin(df.columns)]['IID']
df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)

df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], index=False, sep= '\t', compression= 'gzip')

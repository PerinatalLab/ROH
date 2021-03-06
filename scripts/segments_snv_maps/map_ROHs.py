import pandas as pd
import numpy as np

d= pd.read_csv(snakemake.input[0], delim_whitespace= True)

x= [line.strip() for line in open(snakemake.input[1], 'r')]
fam= [i for i in x if 'fam' in i]
bim= [i for i in x if 'bim' in i]

x= pd.read_csv("".join(bim), delim_whitespace= True, header= None)
fam= pd.read_csv("".join(fam), delim_whitespace=True, header= None)

x.columns= ['CHR','ID', 'CM', 'BP', 'A1','A2']
d= d.loc[:, ['IID','CHR','POS1','POS2']]
x= x.loc[:, ['CHR','BP']]

fam.columns= ['FID', 'IID', 'X1', 'X2', 'X3', 'X4']

CHR= snakemake.wildcards.CHR

z= x.loc[x['CHR']== int(float(CHR)), :]

df_list= list()
chr_df= d.loc[d.CHR== int(float(CHR)),:]
for id in set(chr_df.IID):
	temp_df= chr_df.loc[chr_df.IID== id, :]
	a= z.BP.values
	bh = temp_df.POS2.values
	bl = temp_df.POS1.values
	i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
	temp_df= pd.DataFrame(np.column_stack([z.values[i], temp_df.values[j]]), columns=z.columns.append(temp_df.columns))
#	temp_df= pd.merge(temp_df, x, on= 'CHR', how= 'inner')
#	temp_df= temp_df.loc[(temp_df.BP >= temp_df.POS1) & (temp_df.BP <= temp_df.POS2), :]
	temp_df['Value']= 1
	df_list.append(temp_df)


df= pd.concat(df_list)
df= df.pivot(index= 'BP', columns= 'IID', values= 'Value')
df['BP']= df.index
df['CHR']= CHR


cols = ['CHR', 'BP']  + [col for col in df if col not in 'CHR,BP']
df= df[cols]

o_ids= fam[~fam.IID.isin(df.columns)]['IID']
df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)
df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], index=False, sep= '\t', compression= 'gzip')


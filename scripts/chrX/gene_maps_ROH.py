import pandas as pd
import numpy as np
from functools import reduce

def getOverlap(start0, end0, start1, end1):
        return (np.maximum(0, np.minimum(end0, end1) - np.maximum(start0, start1))) / (end1 - start1)

g= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
d= pd.read_csv(snakemake.input[1], delim_whitespace= True, header= 0)

x= [line.strip() for line in open(snakemake.input[2], 'r')]
fam= [i for i in x if 'fam' in i]

fam= pd.read_csv("".join(fam), delim_whitespace=True, header= None)
fam.columns= ['FID', 'IID', 'X1', 'X2', 'X3', 'X4']

d= d.loc[:, ['IID','CHR','POS1','POS2']]
d.columns= ['IID', 'chr', 'POS1', 'POS2']

df_list= list()

for id in set(d.IID):
	temp_df= d.loc[d.IID== id, :]
	df= pd.merge(temp_df, g, on= ['chr'], how= 'left')
	df['gene']= df.gene.str.replace('.0', '')
	df['overlap']= getOverlap(df.POS1, df.POS2, g.pos1, g.pos2)
	df= df.loc[df.overlap>0, :]
	df['overlap']= 1
	df_list.append(df)

df= pd.concat(df_list)
df.drop_duplicates(['gene', 'IID'], inplace= True)
df= df.pivot(index= 'gene', columns= 'IID', values= 'overlap')
df['gene']= df.index

cols = list(df.columns)
cols = [cols[-1]] + cols[:-1]
df = df[cols]

o_ids= fam[~fam.IID.isin(df.columns)]['IID']
df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)
df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], sep= '\t', index= False, header= True)


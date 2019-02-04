import pandas as pd
import numpy as np
from functools import reduce

def getOverlap(start0, end0, start1, end1):
        return np.maximum(0, np.minimum(end0, end1) - np.maximum(start0, start1))

g= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
d= pd.read_csv(snakemake.input[1], delim_whitespace= True, header= 0)
fam= pd.read_csv(snakemake.input[2], delim_whitespace=True, header= None)
d= d.loc[:, ['IID','CHR','POS1','POS2']]

fam.columns= ['FID', 'IID', 'X1', 'X2', 'X3', 'X4']

if 'rott' in snakemake.input[1]:
        CHR= snakemake.wildcards[1]
if 'harvest' in snakemake.input[1]:
        CHR= snakemake.wildcards[2]


df_list= list()
chr_df= d.loc[d.CHR== int(float(CHR)),:]

for id in set(chr_df.IID):
        temp_df= chr_df.loc[chr_df.IID== id, :]
        temp_df= pd.merge(temp_df, g, left_on= ['CHR'], right_on= ['chr'], how= 'inner')
        temp_df['over']= getOverlap(temp_df.POS1, temp_df.POS2, temp_df.start, temp_df.end)
        temp_df= temp_df.loc[temp_df.over>0, :]
        indiv_map= temp_df.groupby(['chr', 'gene'])['over'].sum()
	indiv_map.rename(columns= {'over': id}, inplace= True)
        indiv_map.set_index(['chr','gene'], inplace= True)
        df_list.append(indiv_map)
df= reduce(lambda df1, df2: df1.merge(df2, "outer", left_index= True, right_index= True), df_list)
df.reset_index(inplace= True)
o_ids= fam[~fam.IID.isin(df.columns)]['IID']
df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)
df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], sep= '\t', index= False, header= True)


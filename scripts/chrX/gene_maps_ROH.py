import pandas as pd
import numpy as np
from functools import reduce

def getOverlap(start0, end0, start1, end1):
        return (np.maximum(0, np.minimum(end0, end1) - np.maximum(start0, start1))) / (end1 - start1)

g= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
g= g.loc[(g.pos2 - g.pos1)>0, :]
g= g.loc[g.chr== 23, :]
temp_list= list()

input_files= [infile for infile in snakemake.input if infile.endswith('.hom')]

for infile in input_files:
        temp_chr= pd.read_csv(infile, delim_whitespace= True, header= 0)
        temp_df= temp_chr.loc[temp_chr.CHR== 23, :]
        temp_list.append(temp_df)

chr_df= pd.concat(temp_list)

df_list= list()

temp_chr= chr_df.loc[chr_df.CHR== 23, :]
for iid in set(chr_df.IID):
        temp_df= temp_chr.loc[temp_chr.IID== iid, :]
        df= pd.merge(temp_df, g, left_on= 'CHR', right_on= 'chr', how= 'inner')
        df['gene']= df.gene.str.replace('.0', '', regex= False)
        df['overlap']= getOverlap(df.POS1, df.POS2, df.pos1, df.pos2)
        df= df.loc[df.overlap>0, :]
        df_list.append(df)

df= pd.concat(df_list)
df= df.groupby(['gene', 'IID'])['overlap'].sum().reset_index()
df['overlap']= np.where(df.overlap > 0, 1, 0)
df= df.pivot(index= 'gene', columns= 'IID', values= 'overlap')
df['gene']= df.index

cols = list(df.columns)
cols = [cols[-1]] + cols[:-1]
df = df[cols]

input_files= [infile for infile in snakemake.input if infile.endswith('.hom.indiv')]
for infile in input_files:
        temp_df= pd.read_csv(infile, delim_whitespace= True, header= 0)
        temp_list.append(temp_df)

iid_df= pd.concat(temp_list)
o_ids= iid_df.loc[~iid_df.IID.isin(df.columns)]['IID']

df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)
df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], sep= '\t', index= False, header= True)



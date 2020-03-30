import pandas as pd
import numpy as np
from functools import reduce



def getOverlap(start0, end0, start1, end1):
        return (np.maximum(0, np.minimum(end0, end1) - np.maximum(start0, start1))) / (end1 - start1)

g= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
g= g.loc[(g.pos2 - g.pos1)>0, :]
lof= pd.read_csv(snakemake.input[1], sep= '\t', header= 0)
temp_list= list()

input_files= [infile for infile in snakemake.input if infile.endswith('.hom')]

def merge_range(df_range, df_val):
	a = df_val.cM.values
	bh = df_range.POS2.values
	bl = df_range.POS1.values
	i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
	stacked_df= pd.DataFrame(np.column_stack([df_val.values[i], df_range.values[j]]), columns= df_val.columns.append(df_range.columns))
	return stacked_df.drop_duplicates(subset= ['IID', 'CHR', 'POS1', 'POS2'], keep= 'first')

for infile in input_files:
	temp_chr= pd.read_csv(infile, delim_whitespace= True, header= 0)
	for CHR in set(temp_chr.CHR):
		temp_df= temp_chr.loc[temp_chr.CHR== CHR, :]
		lof_chr= lof.loc[lof.chr== CHR, :]
		temp_list.append(merge_range(temp_df, lof_chr))

chr_df= pd.concat(temp_list)

df_list= list()

for CHR in set(chr_df.chr):
	temp_chr= chr_df.loc[chr_df.chr== CHR, :]
	for iid in set(chr_df.IID):
		temp_df= temp_chr.loc[temp_chr.IID== iid, :]
		df= pd.merge(temp_df, g, left_on= 'CHR', right_on= 'chr', how= 'inner')
		df['gene']= df.gene.str.replace('.0', '', regex= False)
		df['overlap']= getOverlap(df.POS1, df.POS2, df.pos1, df.pos2)
		df= df.loc[df.overlap>0, :]
	#	df['overlap']= 1
		df_list.append(df)

df= pd.concat(df_list)
#df.drop_duplicates(['gene', 'IID'], inplace= True)
df= df.groupby(['gene', 'IID'])['overlap'].sum().reset_index()
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


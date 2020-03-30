import pandas as pd
import numpy as np

CHR= snakemake.wildcards.CHR

temp_list= list()
input_files= [infile for infile in snakemake.input if infile.endswith('.hom')]
for infile in input_files:
	temp_df= pd.read_csv(infile, delim_whitespace= True, header= 0)
	temp_df= temp_df.loc[temp_df['CHR']== int(float(CHR)),:]
	temp_list.append(temp_df)

chr_df= pd.concat(temp_list)

a= pd.concat([chr_df.POS1, chr_df.POS2])
a= np.unique(a)
a= np.sort(a)

df_list= list()

bh= chr_df.POS2.values
bl= chr_df.POS1.values
i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
df= pd.DataFrame(np.column_stack([a[i], chr_df.values[j]]), columns= ['pos'].append(chr_df.columns))
df.columns= ['start'] + chr_df.columns.values.tolist()

df.sort_values(['IID', 'start'], inplace= True)
df['end']= df.groupby(['IID', 'POS1', 'POS2'])['start'].shift(-1)
df= df[['IID', 'start', 'end']]
df.dropna(inplace= True)
df['Value']= 1
df['segment']= df['start'].astype(str) + ':' + df['end'].astype(str)
df= df.pivot(index= 'segment', columns= 'IID', values= 'Value')
df['segment']= df.index
df['CHR']= CHR

cols = ['CHR', 'segment']  + [col for col in df if col not in 'CHR,segment']
df= df[cols]

input_files= [infile for infile in snakemake.input if infile.endswith('.hom.indiv')]
for infile in input_files:
        temp_df= pd.read_csv(infile, delim_whitespace= True, header= 0)
        temp_list.append(temp_df)

iid_df= pd.concat(temp_list)

o_ids= iid_df.loc[~iid_df.IID.isin(df.columns)]['IID']

df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)

df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], index=False, sep= '\t', compression= 'gzip')

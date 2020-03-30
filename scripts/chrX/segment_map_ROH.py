import pandas as pd
import numpy as np

temp_list= list()
input_files= [infile for infile in snakemake.input[1:] if infile.endswith('.hom')]
for infile in input_files:
        temp_df= pd.read_csv(infile, delim_whitespace= True, header= 0)
        temp_list.append(temp_df)

chr_df= pd.concat(temp_list)
a= pd.concat([chr_df.POS1, chr_df.POS2])
a= np.unique(a)
a= np.sort(a)

df_list= list()

for id in set(chr_df.IID):
	temp_df= chr_df.loc[chr_df.IID== id, :]	
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


input_files= [infile for infile in snakemake.input if infile.endswith('.hom.indiv')]
for infile in input_files:
        temp_df= pd.read_csv(infile, delim_whitespace= True, header= 0)
        temp_list.append(temp_df)

iid_df= pd.concat(temp_list)
o_ids= iid_df.loc[~iid_df.IID.isin(df.columns)]['IID']

df= pd.concat((df, pd.DataFrame(columns= o_ids, index= df.index)), axis= 1)

df.fillna(0, inplace= True)

df.to_csv(snakemake.output[0], index=False, sep= '\t', compression= 'gzip')

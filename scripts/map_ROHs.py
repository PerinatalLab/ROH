import pandas as pd

d= pd.read_csv(snakemake.input[0], delim_whitespace= True)
x= pd.read_csv(snakemake.input[1], delim_whitespace= True, header= None)
fam= pd.read_csv(snakemake.input[2], delim_whitespace=True, header= None)

x.columns= ['CHR','ID', 'CM', 'BP', 'A1','A2']
d= d.loc[:, ['FID','CHR','POS1','POS2']]
x= x.loc[:, ['CHR','BP']]

fam.columns= ['FID', 'IID', 'X1', 'X2', 'X3', 'X4']

CHR= snakemake.wildcards[2]
df_list= list()
chr_df= d.loc[d.CHR== int(float(CHR)),:]
for id in set(chr_df.FID):
	temp_df= chr_df.loc[chr_df.FID== id,:]
	temp_df= pd.merge(temp_df, x, on= 'CHR') 
	temp_df= temp_df.loc[(temp_df.BP >= temp_df.POS1) & (temp_df.BP <= temp_df.POS2), :]
	temp_df[id]= 1
	temp_df.set_index(['CHR', 'BP'], inplace=True)
	temp_df.drop(['POS1', 'POS2', 'FID'], inplace= True, axis= 1)
	df_list.append(temp_df)
df= pd.concat(df_list, axis=1)
o_ids= fam[~fam.FID.isin(df.columns)]['FID']
df= pd.concat((df, pd.DataFrame(columns= o_ids, index= temp_df.index)), axis= 1)
df.fillna(0, inplace= True)
df= df.reset_index(level= [0,1])
df.to_csv(snakemake.output[0], index=False, sep= '\t')


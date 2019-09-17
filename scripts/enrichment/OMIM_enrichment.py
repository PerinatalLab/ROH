import numpy as np
import pandas as pd
import scipy.stats as st

def getOverlap(start0, end0, start1, end1):
	return (np.maximum(0, np.minimum(end0, end1) - np.maximum(start0, start1))) / (np.array(end0) - np.array(start0) )

x= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)

omim= pd.read_csv(snakemake.input[1], sep= '\t', header= 0)


df_list= list()
for CHR in set(x.chr_NC):
	temp_df= x.loc[x.chr_NC== CHR, :]
	temp_df['newid']= range(temp_df.shape[0])
	temp_df2= pd.merge(temp_df, omim, left_on= ['chr_NC'], right_on= ['chr'], how= 'inner')
	temp_df2['overlap_NC']= getOverlap(temp_df2.cM1_NC, temp_df2.cM2_NC, temp_df2.cM1_OMIM, temp_df2.cM2_OMIM)
	newdf= temp_df2.groupby(temp_df2.newid).agg({'overlap_NC': np.sum})
	newdf.reset_index(inplace= True)
	temp_df= pd.merge(newdf, temp_df, on= ['newid'])
	df_list.append(temp_df)

df_NC= pd.concat(df_list)
df_NC['overlap_NC']= np.where(df_NC.overlap_NC> 1, 1, df_NC.overlap_NC)
df_NC= df_NC.loc[:, ['segment_NC', 'segment', 'overlap_NC']]

df_list= list()
for CHR in set(x.chr_NC):
	temp_df= x.loc[x.chr_NC== CHR, :]
	temp_df= temp_df.groupby('segment').head(1)
	temp_df['newid']= range(temp_df.shape[0])
	temp_df2= pd.merge(temp_df, omim, on= ['chr'], how= 'inner')
	temp_df2['overlap_HC']= getOverlap(temp_df2.cM1, temp_df2.cM2, temp_df2.cM1_OMIM, temp_df2.cM2_OMIM)
	newdf= temp_df2.groupby(temp_df2.newid).agg({'overlap_HC': np.sum})
	newdf.reset_index(inplace= True)
	temp_df= pd.merge(newdf, temp_df, on= ['newid'])
	df_list.append(temp_df)

df_HC= pd.concat(df_list)
df_HC['overlap_HC']= np.where(df_HC.overlap_HC> 1, 1, df_HC.overlap_HC)
df_HC= df_HC.groupby('segment').head(1)
df_HC.drop(['segment_NC'], inplace= True, axis= 1)

df= pd.merge(df_HC, df_NC, on= 'segment')
df['bigger_overlap']= np.where(df.overlap_HC>= df.overlap_NC, 1, 0)
df.drop_duplicates(subset= ['segment', 'segment_NC'], keep= 'first', inplace= True)
newdf= pd.DataFrame(df.groupby('segment')['bigger_overlap'].mean().reset_index())
newdf['bigger_overlap']= 1 - newdf['bigger_overlap']
newdf.to_csv(snakemake.output[0], sep= '\t', header= True, index= False)

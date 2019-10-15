#!/usr/bin/python3

import pandas as pd
import numpy as np

def pheno():
        d= pd.read_csv(snakemake.input[0], delim_whitespace= True, header= 0)
        dX= pd.read_csv(snakemake.input[1], delim_whitespace= True, header= 0)
	dX= d[['IID', 'KB', 'KBAVG', 'NSEG']]
	dX.columns= ['IID', 'KBX', 'KBAVGX', 'NSEGX']
        bim= [line.strip() for line in open(snakemake.input[2], 'r')]
        bim= "".join(bim[1])
        bim= pd.read_csv(bim, sep= '\t', header= None, names= ['chr', 'snp', 'cM', 'pos', 'A1', 'A2'])

        bp= bim.groupby(['chr'])['pos'].diff(1).sum() / 1000000

        d= pd.merge(d, dX, on= ['IID'], how= 'outer')
	d['KB']= d[['KB', 'KBX']].sum(axis= 1)
	d['KBAVG']= d[['KBAVG', 'KBAVGX']].sum(axis= 1)
	d['NSEG']= d[['NSEG', 'NSEGX']].sum(axis= 1)
        d['KB']= d['KB'] / 1000000 * 1000
        d['FKB']= d['KB'].divide(bp, axis=1)
        return d

d= pheno()
d.to_csv(snakemake.output[0], sep= '\t', index= False)

import pandas as pd
import numpy as np
import os
import gzip
from functools import reduce
import scipy.stats as st
import statsmodels.stats.multitest as multi

cohort_nms= ['harvestm12', 'harvestm24','rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay']
smpl_nms= ['maternal','paternal', 'fetal']
batch_nms= ['m12', 'm24']
CHR_nms= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22]

# Other arguments:

pruning_nms= ['none', 'soft', 'moderate', 'hard']

dens_nms= [5]
SNP_nms= [15, 25, 50, 75, 100, 150, 200, 350, 500]
length_nms= [0.0000001]
het_nms= [0, 1]
GAP_nms= [5]

dens_bp= [5000]
SNP_bp= [15, 25, 50, 75, 100, 150, 200, 350, 500]
length_bp= [0.0000001]
het_bp= [0, 1]
GAP_bp= [5000]

# Functions

def isfloat(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

# Rules

rule all:
        'Collect the main outputs of the workflow.'
        input:
                expand('/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt', cohort= cohort_nms, sample= smpl_nms),
                expand('/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',cohort= cohort_nms),
                expand('/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt', cohort= cohort_nms),
#		expand('reports/ROH_{cohort}_analysis.html', cohort= cohort_nms),
                'figures/figure2.eps',
                'figures/figure1.eps',
                'figures/S1_figure.eps',
                'figures/S2_figure.eps',
		'tables/S1_table.txt',
                'figures/figureX.eps'

include: 'scripts/cox/Snakefile'
include: 'scripts/figures/Snakefile'
include: 'scripts/metaanalysis/Snakefile'
include: 'scripts/segments_snv_maps/Snakefile'
include: 'scripts/enrichment/Snakefile'
include: 'scripts/frequency/Snakefile'
include: 'scripts/reports/Snakefile'
include: 'scripts/imputed/Snakefile'
include: 'scripts/phasing/Snakefile'
include: 'scripts/ROH_calling/Snakefile'

## Snakemake code

rule ids_to_keep:
	'List maternal, paternal and fetal ids acceptable by PLINK for --keep.'
	input:
		'/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv'
	output:
                '/mnt/work/pol/ROH/{cohort}/pheno/maternal_ids',
                '/mnt/work/pol/ROH/{cohort}/pheno/paternal_ids',
                '/mnt/work/pol/ROH/{cohort}/pheno/fetal_ids',
		'/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt'
	run:
		if 'harvest' in wildcards.cohort:
			d= pd.read_csv(input[0], sep= '\t', header= 0)
			d.dropna(subset= ['Role'], inplace= True)
			x= d.pivot(index='PREG_ID_1724', columns='Role', values= [ 'SentrixID_1'])
			x.columns= x.columns.droplevel()
			x.reset_index(inplace=True)
			x.columns= ['PREG_ID_1724', 'Child', 'Father', 'Mother']
			x.dropna(inplace= True)
			x.to_csv(output[0], header= None, columns= ['Mother', 'Mother'], index= False, sep= '\t')
	                x.to_csv(output[2], header= None, columns= ['Child', 'Child'], index= False, sep= '\t')
		        x.to_csv(output[1], header= None, columns= ['Father', 'Father'], index= False, sep= '\t')
			x.to_csv(output[3], header= True, sep= '\t', index= False)
		if 'harvest' not in wildcards.cohort:
			d= pd.read_csv(input[0], delim_whitespace= True, header= 0)
			d.dropna(subset= ['Role'], inplace= True)
			x= d.pivot(index= 'PREG_ID_315', columns= 'Role', values= ['FID', 'SentrixID'])
			x.columns= x.columns.droplevel()
			x.iloc[:,2]= x.iloc[:,2].fillna(x.iloc[:,0])
			x.iloc[:,2]= x.iloc[:,2].fillna(x.iloc[:,1])
			x= x.iloc[:, 2:]
			x.reset_index(inplace=True)
			x.columns= ['PREG_ID_315', 'FID', 'Child', 'Father', 'Mother']
			x['PREG_ID_315']= x['PREG_ID_315'].astype(int)
			x= x[pd.to_numeric(x['FID'], errors='coerce').notnull()]
			x['FID']= x['FID'].astype(int)
			x.dropna(inplace= True)
			x.to_csv(output[0], header= None, columns= ['FID', 'Mother'], index= False, sep= '\t')
		        x.to_csv(output[2], header= None, columns= ['FID', 'Child'], index= False, sep= '\t')
			x.to_csv(output[1], header= None, columns= ['FID', 'Father'], index= False, sep= '\t')
			x.to_csv(output[3], header= True, sep= '\t', index= False)

rule phenofile:
        'Merge all data necessary to create a phenotype file with ROH.'
        input:
                '/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom',
                '/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom.indiv',
                '/mnt/work/pol/{cohort}/pheno/{cohort}_mfr.csv',
                '/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv',
                '/mnt/work/pol/{cohort}/pca/{cohort}_pca.txt',
                '/mnt/work/pol/{cohort}/relatedness/all_{cohort}.kin0',
		'/mnt/archive/HARVEST/delivery-fhi/data/genotyped/m12/m12-genotyped.fam',
		'/mnt/work/pol/ROH/{cohort}/runs/{sample}_input_ROH_geno.txt',
		'/mnt/work/pol/{cohort}/pheno/flag_list.txt',
		'/mnt/work/pol/{cohort}/pca/all_pca_outliers_hapmap.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{{sample}}.bim', pruning= pruning_nms)
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
        script:
                'scripts/pheno_file.py'

rule dl_genetic_map:
	'Download the genetic map estimated in 1KG (https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html), from IMPUTE2.'
	output:
		temp(expand('/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_chr{CHR}_combined_b37.txt', CHR= CHR_nms))
	shell:
		'''
		wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -P /mnt/work/pol/ROH/1KG/
		tar -xvzf /mnt/work/pol/ROH/1KG/1000GP_Phase3.tgz
		mv 1000GP_Phase3 /mnt/work/pol/ROH/1KG/
		rm /mnt/work/pol/ROH/1KG/1000GP_Phase3.tgz /mnt/work/pol/ROH/1KG/1000GP_Phase3/*hap.gz /mnt/work/pol/ROH/1KG/1000GP_Phase3/*.legend.gz /mnt/work/pol/ROH/1KG/1000GP_Phase3/1000GP_Phase3.sample
		'''

rule concat_genetic_map:
	'Concat all genetic map files.'
	input:
		expand('/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_chr{CHR}_combined_b37.txt', CHR= CHR_nms)
	output:
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt'
	run:
		df= pd.DataFrame()
		for infile in input:
			d= pd.read_csv(infile, header= 0, sep= ' ')
			x= infile.split('_')
			x= [s.replace('chr', '') for s in x if s.replace('chr','').isdigit()]
			x= ''.join(x)
			d['chr']= x
			d= d[['chr', 'position', 'COMBINED_rate(cM/Mb)', 'Genetic_Map(cM)']]
			df= df.append(d)
		df.to_csv(output[0], sep= ' ', index= False, header= True)

rule map_format_geneticmap:
	'Format the genetic map into .map file format from PLINK'
	input:
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt',
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.haps'
	output:
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/{cohort}_genetic_map_combined_b37_{CHR}.map'
	run:
		d= pd.read_csv(input[0], sep= ' ', header= 0)
		d['SNP']= d.chr.map(str) + ':' + d.position.map(str)
		d= d[['chr', 'SNP', 'Genetic_Map(cM)', 'position']]
		d= d.loc[d.chr== wildcards.CHR, :]
		ilu= np.loadtxt(input[1], usecols= (0,2), delimiter= ' ')
		ilu= pd.DataFrame({'chr': ilu[:,0], 'position': ilu[:,1]})
		df= pd.merge(d, ilu, on= ['chr', 'position'], how= 'right')
		df= df.loc[df['Genetic_Map(cM)'].isna(), :]	
		df['Genetic_Map(cM)']= np.interp(df.position, d['position'], d['Genetic_Map(cM)'])
		d= d.append(df)
		d['SNP']= d.chr.map(str) + ':' + d.position.map(str)
		d.to_csv(output[0], header= False, index= False, sep= '\t')


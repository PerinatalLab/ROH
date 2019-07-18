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

include: 'scripts/Snakefile'

rule all:
	'Collect the main outputs of the workflow.'
	input:
		expand('/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt', cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt', cohort= cohort_nms),	
		expand('reports/ROH_{cohort}_analysis.html', cohort= cohort_nms)

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

rule copy_harvest_genotyped:
	''
	input:
		expand('/mnt/archive/HARVEST/delivery-fhi/data/genotyped/m12/m12-genotyped.{ext}', ext= ['bed', 'bim', 'fam']),
                expand('/mnt/archive/HARVEST/delivery-fhi/data/genotyped/m24/m24-genotyped.{ext}', ext= ['bed', 'bim', 'fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/harvestm12/genotypes/temp/harvestm12_genotyped.{ext}', ext= ['bed', 'bim', 'fam'])),
                temp(expand('/mnt/work/pol/ROH/harvestm24/genotypes/temp/harvestm24_genotyped.{ext}', ext= ['bed', 'bim', 'fam']))
	shell:
		'''
		cp {input[0]} {output[0]}; cp {input[1]} {output[1]}; cp {input[2]} {output[2]}
                cp {input[3]} {output[3]}; cp {input[4]} {output[4]}; cp {input[5]} {output[5]}
		'''
rule copy_rott_genotyped:
        'Copy genotype PLINK files and change name.'
        input:
                expand('/mnt/archive/ROTTERDAM1/delivery-fhi/data/genotyped/genotyped.{ext}', ext= ['bed', 'bim', 'fam']),
                expand('/mnt/archive/ROTTERDAM2/delivery-fhi/data/genotyped/genotyped.{ext}', ext= ['bed', 'bim', 'fam'])
        output:
                temp(expand('/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_genotyped.{ext}', ext= ['bed', 'bim', 'fam'])),
                temp(expand('/mnt/work/pol/ROH/rotterdam2/genotypes/temp/rotterdam2_genotyped.{ext}', ext= ['bed', 'bim', 'fam'])),
        shell:
                '''
                cp {input[0]} {output[0]}; cp {input[1]} {output[1]}; cp {input[2]} {output[2]}
                cp {input[3]} {output[3]}; cp {input[4]} {output[4]}; cp {input[5]} {output[5]}
                '''

rule copy_norment_genotyped:
        'Copy genotype PLINK files and change name.'
        input:
                expand('/mnt/archive/NORMENT1/delivery-fhi/data/genotyped/feb18/genotyped.{ext}', ext= ['bed', 'bim', 'fam']),
                expand('/mnt/archive/NORMENT1/delivery-fhi/data/genotyped/may16/genotyped.{ext}', ext= ['bed', 'bim', 'fam'])
        output:
                temp(expand('/mnt/work/pol/ROH/normentfeb/genotypes/temp/normentfeb_genotyped.{ext}', ext= ['bed', 'bim', 'fam'])),
                temp(expand('/mnt/work/pol/ROH/normentmay/genotypes/temp/normentmay_genotyped.{ext}', ext= ['bed', 'bim', 'fam']))
        shell:
                '''
                cp {input[0]} {output[0]}; cp {input[1]} {output[1]}; cp {input[2]} {output[2]}
                cp {input[3]} {output[3]}; cp {input[4]} {output[4]}; cp {input[5]} {output[5]}
                '''


rule exclude_non_biallelic:
        'Set range file for multi-allelic SNPs.'
        input:
                '/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_genotyped.bim'
        output:
                temp('/mnt/work/pol/ROH/{cohort}/multiallelic.txt')
        run:
                d= pd.read_csv(input[0], sep= '\t', header= None)
                d.columns= ['chr', 'SNP', 'X', 'pos', 'A1', 'A2']
		d= d[d[['chr', 'pos']].duplicated(keep=False)]
		d= d[['chr', 'pos', 'pos', 'X']]
		d.to_csv(output[0], sep= '\t', header= False, index= False)

rule split_bed:
        'Modify the bed file: remove CHR 23, 24, 25 and 26, maf <=0.05 and split file by sample.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_genotyped.{ext}', ext= ['bed','bim','fam']),
                '/mnt/work/pol/ROH/{cohort}/pheno/{sample}_ids',
                '/mnt/work/pol/ROH/{cohort}/multiallelic.txt'
        output:
                temp(expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'log']))
        params:
                '/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_genotyped',
                '/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_genotyped_{sample}'
        shell:
                '~/soft/plink --bfile {params[0]} --exclude range {input[4]} --maf 0.05 --keep {input[3]} --make-bed --chr 1-22 --make-founders --out {params[1]}'

rule multi_pruning:
	'Filter PLINK file according to different pruning parameters.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'log'])	
	output:
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/soft/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['prune.out','prune.in', 'log'])),
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/moderate/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['prune.out','prune.in', 'log'])),
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/hard/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['prune.out','prune.in', 'log']))
	params:
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_genotyped_{sample}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/soft/{cohort}_genotyped_{sample}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/moderate/{cohort}_genotyped_{sample}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/hard/{cohort}_genotyped_{sample}'
	shell:
		"""
		~/soft/plink --bfile {params[0]} --indep-pairwise 50 5 0.9 --out {params[1]}
		~/soft/plink --bfile {params[0]} --indep-pairwise 50 5 0.5 --out {params[2]}
		~/soft/plink --bfile {params[0]} --indep-pairwise 50 5 0.1 --out {params[3]}
		"""

rule move_none_pruning:
	'Move PLINK files not pruned to wildcard.pruning == none folder.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['bed','bim','fam'])
	output:
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/none/pruned{{cohort}}_{{sample}}.{ext}', ext= ['bed','bim','fam'])
	params:
		'/mnt/work/pol/ROH/{cohort}/genotypes/none/'
	shell:
		"""
		mkdir -p {params[0]}
		cp {input[0]} {output[0]}
		cp {input[1]} {output[1]}
		cp {input[2]} {output[2]}
		"""

rule plink_bfile_prune:
	'Exclude genetic variants in prune.out files (obtained with rule plink_split_bed).'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/hard/{cohort}_genotyped_{sample}.prune.out',
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/soft/{cohort}_genotyped_{sample}.prune.out',
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/moderate/{cohort}_genotyped_{sample}.prune.out',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_genotyped_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam'])
	output:
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/hard/pruned{{cohort}}_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'log']),
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/soft/pruned{{cohort}}_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'log']),
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/moderate/pruned{{cohort}}_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'log'])
	params:
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_genotyped_{sample}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/hard/pruned{cohort}_{sample}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/soft/pruned{cohort}_{sample}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/moderate/pruned{cohort}_{sample}'
	shell:
		'''
		~/soft/plink --bfile {params[0]} --exclude {input[0]} --make-bed --out {params[1]}
		~/soft/plink --bfile {params[0]} --exclude {input[1]} --make-bed --out {params[2]}
		~/soft/plink --bfile {params[0]} --exclude {input[2]} --make-bed --out {params[3]}
		'''

rule replace_bp_cm:
	'PLINK cannot use cM to estimate ROH length, so we replace bp position to cM in the .bim file.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_{sample}.bim',
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt'
	output:
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/cm_pruned{cohort}_{sample}.bim'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= None)
		d.columns= ['chr', 'SNP', 'X', 'pos', 'A1', 'A2']
		g= pd.read_csv(input[1], sep= ' ', header= 0)
		g= g[['chr', 'Genetic_Map(cM)', 'position']]
		g.columns= ['chr', 'Genetic_Map(cM)', 'pos']
		df= pd.merge(d, g, on= ['chr', 'pos'], how= 'left')
		df_miss= df.loc[df['Genetic_Map(cM)'].isna(), :]
		newdf= pd.DataFrame()
		for CHR in set(df_miss.chr):
			df_temp= df_miss.loc[df_miss.chr== CHR, :]
			g_temp= g.loc[g.chr== CHR, :]
			df_temp['newX']= np.interp(df_temp['pos'], g_temp['pos'], g_temp['Genetic_Map(cM)'])
			newdf= newdf.append(df_temp)
		newdf= newdf[['chr','pos', 'newX']]
		df= pd.merge(df, newdf, on= ['chr', 'pos'], how= 'left')
		df['X']= np.where(df['Genetic_Map(cM)'].isna(), df['newX'], df['Genetic_Map(cM)'])
		df['X']= (df.X*10**4).round() * 100
		df['X']= df['X'] + df.groupby(['chr', 'X']).cumcount()
		df[['pos', 'X']]= df[['X','pos']]
		df= df[['chr', 'SNP', 'X', 'pos', 'A1', 'A2']]
		df.to_csv(output[0], sep= '\t', header= False, index= False)

rule run_ROH_multi_arg:
	'Estimate ROH using multiple arguments.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_fetal.bed',
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/cm_pruned{cohort}_fetal.bim',
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_fetal.fam'
	output:
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/multi/{{pruning}}_fetal_{{dens}}_{{SNP}}_{{length}}_{{het}}_{{GAP}}.{ext}', ext= ['log', 'hom', 'hom.summary'])),
		'/mnt/work/pol/ROH/{cohort}/multi/{pruning}_fetal_{dens}_{SNP}_{length}_{het}_{GAP}.hom.indiv'
	params:
		'/mnt/work/pol/ROH/{cohort}/multi/{pruning}_fetal_{dens}_{SNP}_{length}_{het}_{GAP}'
	run:
		SNPwm= round(float(wildcards.SNP) * 0.05)
		GAP= int(float(wildcards.GAP) * 1000)
		dens= int(float(wildcards.dens) * 1000)
		shell("/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --homozyg-window-snp {wildcards.SNP} --homozyg-snp {wildcards.SNP} --homozyg-kb 0.0000001 --homozyg-gap {GAP} --homozyg-window-missing {SNPwm} --homozyg-window-threshold 0.0005 --homozyg-window-het {wildcards.het} --homozyg-density {dens} --out {params[0]}")

rule run_ROH_multi_arg_bp:
        'Estimate ROH using multiple arguments.'
        input:
                '/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_fetal.bed',
                '/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_fetal.bim',
                '/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_fetal.fam'
        output:
                temp(expand('/mnt/work/pol/ROH/{{cohort}}/multi/{{pruning}}_bpfetal_{{densbp}}_{{SNPbp}}_{{lengthbp}}_{{hetbp}}_{{GAPbp}}.{ext}', ext= ['log', 'hom', 'hom.summary'])),
		'/mnt/work/pol/ROH/{cohort}/multi/{pruning}_bpfetal_{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}.hom.indiv'
        params:
                '/mnt/work/pol/ROH/{cohort}/multi/{pruning}_bpfetal_{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}'
	run:
		SNPwm= round(float(wildcards.SNPbp) * 0.05)
		shell("/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --homozyg-window-snp {wildcards.SNPbp} --homozyg-snp {wildcards.SNPbp} --homozyg-kb {wildcards.lengthbp} --homozyg-gap {wildcards.GAPbp} --homozyg-window-missing {SNPwm} --homozyg-window-threshold 0.0005 --homozyg-window-het {wildcards.hetbp} --homozyg-density {wildcards.densbp} --out {params[0]}")


                '/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_pca.txt',
                '/mnt/work/pol/ROH/{cohort}/ibd/to_phase.fam',
                '/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_fetal.fam',
                '/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt',
                '/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_fetal.txt'


rule determine_arguments_ROH:
	'Determine ROH estimation arguments that maximise ROH - parental IBD correlation.'
	input:
                '/mnt/work/pol/{cohort}/pca/{cohort}_pca.txt',
                '/mnt/work/pol/ROH/{cohort}/ibd/to_phase.fam',
                '/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_fetal.fam',
                '/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt',
		'/mnt/work/pol/{cohort}/pheno/flag_list.txt',
		'/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt',
		'/mnt/work/pol/{cohort}/relatedness/all_{cohort}.kin0',
		'/mnt/work/pol/{cohort}/pca/all_pca_outliers_hapmap.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/multi/{pruning}_fetal_{dens}_{SNP}_{length}_{het}_{GAP}.hom.indiv', dens= dens_nms, SNP= SNP_nms, length= length_nms, het= het_nms, GAP= GAP_nms, pruning= pruning_nms),	
		expand('/mnt/work/pol/ROH/{{cohort}}/multi/{pruning}_bpfetal_{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}.hom.indiv', densbp= dens_bp, SNPbp= SNP_bp, lengthbp= length_bp, hetbp= het_bp, GAPbp= GAP_bp, pruning= pruning_nms)
	output:
		'/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',
		'/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt'
	script:
		'scripts/R2_ROH_IBD.R'

rule estimate_ROH:
	'''
	Obtain ROH estimates using PLINK 1.9.
	Configuration according to file "/mnt/work/pol/ROH/arguments/max_R2.txt"
	'''
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{{sample}}.{ext}', pruning= pruning_nms, ext= ['bed', 'bim', 'fam']),
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/cm_pruned{{cohort}}_{{sample}}.bim', pruning= pruning_nms),
		'/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt'
	output:
		'/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom.indiv',
		'/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom',
		'/mnt/work/pol/ROH/{cohort}/runs/{sample}_input_ROH_geno.txt',
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/runs/{{cohort}}_{{sample}}.{ext}', ext= ['log', 'hom.summary']))
	params:
		'/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}'
	run:
		maxim= [df for df in input if 'max_R2' in df]
		parlist= [line.strip() for line in open("".join(maxim), 'r')]
		parlist= [float(x) for x in parlist]
		if parlist[0]== 0:
			prun= 'none'
		if parlist[0]== 1:
			prun= 'soft'
		if parlist[0]== 2:
			prun= 'moderate'
		if parlist[0]== 3:
			prun= 'hard'
		bed= [bed for bed in input if prun in bed and 'bed' in bed]
		fam= [fam for fam in input if prun in fam and 'fam' in fam]
		GAP= round(parlist[5] * 1000)
		SNPwm= round(parlist[3] * 0.05)
		dens= round(parlist[1] * 1000)
		if parlist[1] < 100:
			bim= [bim for bim in input if prun in bim and 'bim' in bim and 'cm' in bim]
			shell("/home/pol.sole.navais/soft/plink --bed {bed} --bim {bim} --fam {fam} --homozyg-window-snp {parlist[2]} --homozyg-snp {parlist[2]} --homozyg-kb {parlist[3]} --homozyg-gap {GAP} --homozyg-window-missing {SNPwm} --homozyg-window-threshold 0.0005 --homozyg-window-het {parlist[4]} --homozyg-density {dens} --out {params}")
		if parlist[1] > 100:
			bim= [bim for bim in input if prun in bim and 'bim' in bim and 'cm' not in bim]
			shell("/home/pol.sole.navais/soft/plink --bed {bed} --bim {bim} --fam {fam} --homozyg-window-snp {parlist[2]} --homozyg-snp {parlist[2]} --homozyg-kb {parlist[3]} --homozyg-gap {parlist[5]} --homozyg-window-missing {SNPwm} --homozyg-window-threshold 0.0005 --homozyg-window-het {parlist[4]} --homozyg-density {parlist[1]} --out {params}")
		l= list(["".join(bed), "".join(bim), "".join(fam)])
		with open(output[2], 'w') as f:
			f.writelines( "%s\n" % item for item in l)

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

rule mapping_ROH_segments:
        'Obtain matrix (rows= segment, columns = subject), with all minimum segmental ROHs per subject (1= homozygous part of ROH).'
        input:
                '/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom',
		'/mnt/work/pol/ROH/{cohort}/runs/{sample}_input_ROH_geno.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{{sample}}.fam', pruning= pruning_nms)
        output:
                temp('/mnt/work/pol/ROH/{cohort}/genotypes/maps/{sample}/segments_maps_{sample}_chr{CHR}.txt.gz')
        script:
                'scripts/segment_map_ROH.py'

rule ROH_freq:
        'Count per-position relative frequency of ROHs.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/maps/{{sample}}/segments_maps_{{sample}}_chr{CHR}.txt.gz', CHR= CHR_nms)
        output:
                '/mnt/work/pol/ROH/{cohort}/runs/frequency/ROH_frequency_{sample}'
        run:
                for i in input:
                        for chunk in pd.read_csv(gzip.open(i), sep ='\t', header= 0, chunksize= 500):
                                chunk.fillna(0, inplace= True)
                                x= chunk.iloc[:,2:].mean(axis= 1)
                                x= pd.concat([chunk.iloc[:,0:2], x], axis= 1, ignore_index= True, sort= False)
                                x.to_csv(output[0], mode= 'a', sep= '\t', header= False, index= False)

rule cox_ph_analysis:
	'Cox proportional hazard analyis on bined ROH calls.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/maps/{sample}/segments_maps_{sample}_chr{CHR}.txt.gz',
		'/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}')
	script:
		'scripts/cox_segments.R'

rule concat_cox:
	'Concat cox results from multiple chromosomes.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/results/maps_cox/{{sample}}/cox_spont_{{sample}}_chr{CHR}', CHR= CHR_nms)
	output:
		temp('/mnt/work/pol/ROH/{cohort}/results/maps_cox/{sample}/cox_spont_{sample}')
	shell:
		'cat {input} > {output[0]}'

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

rule excess_homozygosity:
	'Compute excess homozygosity using PLINK --het flag.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_{sample}.bed',
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_{sample}.bim',
		'/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_{sample}.fam'
	output:
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/results/het/{{pruning}}_excess_hom_{{sample}}.{ext}', ext= ['het', 'log']))
	params:
		'/mnt/work/pol/ROH/{cohort}/results/het/{pruning}_excess_hom_{sample}'
	shell:
		'~/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --het --out {params[0]}'

rule merge_homozygosity:
	'Merge excess homozygosity from different pruning.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/results/het/{pruning}_excess_hom_{{sample}}.het', pruning= pruning_nms)
	output:
		'/mnt/work/pol/ROH/{cohort}/results/het/{sample}_excess_hom.txt'
	run:
		dflist= list()
		for i in input:
			d= pd.read_csv(i, delim_whitespace=True, header=0)
			d= d[['IID', 'F']]
			if 'none' in i:	d.columns= ['IID', 'none_F']
			if 'soft' in i: d.columns= ['IID', 'soft_F']
			if 'moderate' in i: d.columns= ['IID', 'moderate_F']
			if 'hard' in i: d.columns= ['IID', 'hard_F']
			dflist.append(d)
		d= reduce(lambda x, y: pd.merge(x, y, on= 'IID'), dflist)
		d.to_csv(output[0], index=False, header= True, sep= '\t')


rule fam_to_phasing:
        'Obtain a fam file for those samples in which IBD detection was applied.'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/to_imputation/m12/m12-ready-for-imputation.fam',
                '/mnt/archive/HARVEST/delivery-fhi/data/to_imputation/m24/m24-ready-for-imputation.fam',
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/to_phasing/merged/hrc-update-complete-all.fam',
                '/mnt/archive/ROTTERDAM2/delivery-fhi/data/to_phasing/merged/hrc-update-complete.fam',
		'/mnt/archive/NORMENT1/delivery-fhi/data/to_phasing/feb18/merged/hrc-update-complete.fam',
		'/mnt/archive/NORMENT1/delivery-fhi/data/to_phasing/may16/merge/hrc-update-complete.fam'
	output:
		'/mnt/work/pol/ROH/{cohort}/ibd/to_phase.fam'
	run:
		if wildcards.cohort== 'harvestm12': fam= input[0]
		if wildcards.cohort== 'harvestm24': fam= input[1]
		if wildcards.cohort== 'rotterdam1': fam= input[2]
		if wildcards.cohort== 'rotterdam2': fam= input[3]
		if wildcards.cohort== 'normentfeb': fam= input[4]
		if wildcards.cohort== 'normentmay': fam= input[5]
		shell('cp {fam} {output[0]}')

rule split_segments:
	'Split overlapping segments between sub-cohorts.'
	input:
		'/mnt/work/pol/ROH/{cohort}/results/maps_cox/{sample}/cox_spont_{sample}',
		expand('/mnt/work/pol/ROH/{cohort}/results/maps_cox/{{sample}}/cox_spont_{{sample}}', cohort= cohort_nms)
	output:
		temp('/mnt/work/pol/ROH/{cohort}/results/cox_spont_{sample}_temp')
	script:
		'scripts/overlap_split_segment.py'

rule cox_cM_to_bp:
        'Convert cM to bp in cox results.'
        input:
                '/mnt/work/pol/ROH/{cohort}/results/cox_spont_{sample}_temp',
                expand('/mnt/work/pol/ROH/{cohort}/runs/{{sample}}_input_ROH_geno.txt', cohort= cohort_nms),
                expand('/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_{{sample}}.bim', cohort= cohort_nms, pruning= pruning_nms)
        output:
                '/mnt/work/pol/ROH/{cohort}/results/cox_spont_{sample}'
	run:
		d= pd.read_csv(input[0], header= 0, sep= '\t')
		flist= [infile for infile in input if 'geno.txt' in infile]
		bim_list= list()
		for infile in flist:
			x= [line.strip() for line in open("".join(infile), 'r')]
			bim= [i for i in x if 'bim' in i]
			bim= pd.read_csv("".join(bim), delim_whitespace= True, header= None, names= ['chr', 'id', 'pos', 'cM', 'A1', 'A2'])
			bim_list.append(bim)
		bim= pd.concat(bim_list)
		bim.drop_duplicates(subset= ['chr', 'pos'], keep= 'first', inplace= True)
		bim.drop_duplicates(subset= ['chr', 'cM'], keep= 'first', inplace= True)
		bim= bim[['chr', 'pos', 'cM']]
		d= pd.merge(d, bim, left_on= ['chr', 'cM1'], right_on= ['chr', 'cM'])
		d.rename(columns= {'pos': 'pos1'}, inplace= True)
		d= pd.merge(d, bim, left_on= ['chr', 'cM2'], right_on= ['chr', 'cM'])
		d.rename(columns= {'pos': 'pos2'}, inplace= True)
		d.drop(['cM_x', 'cM_y'], axis= 1, inplace= True)
		d.to_csv(output[0], sep= '\t', header= True, index= False)

rule high_conf_segments:
	'Obtain a set of high confidence intervals with an FDR< 0.05 and <0.1.'
	input:
		expand('/mnt/work/pol/ROH/{cohort}/results/cox_spont_{{sample}}', cohort= cohort_nms)
	output:
		'/mnt/work/pol/ROH/results/HC_{sample}_cox_spont',
		'/mnt/work/pol/ROH/results/LC_{sample}_cox_spont',
		'/mnt/work/pol/ROH/results/NC_{sample}_cox_spont'
	run:
		i= 0
		df_list= list()
		for infile in input:
			d= pd.read_csv(infile, sep= '\t', header= 0)
			d['w_zscore']= (d.beta / d.sd) * (1/d.sd)
			d['denominator']= (1/ d.sd)**2
			d['segment']= d.chr.map(str) + ':' + d.cM1.map(str) + ':' + d.cM2.map(str)
			d.rename(columns={'w_zscore': 'w_zscore'+str(i), 'denominator': 'denominator'+ str(i), 'pos1': 'pos1'+ str(i), 'pos2': 'pos2'+ str(i)}, inplace=True)
			d.drop(['cM1', 'cM2', 'chr', 'n', 'beta', 'sd', 'pvalue', 'R', 'R_pvalue'], inplace= True, axis= 1)
			i+= 1
			df_list.append(d)
		df= reduce(lambda x, y: pd.merge(x, y, on = 'segment', how= 'outer'), df_list)
		cols= [col for col in df.columns if 'w_zscore' in col]
		df['numerator']= df[cols].sum(axis= 1)
		cols2= [col for col in df.columns if 'denominator' in col]
		df['denominator']= df[cols2].sum(axis =1)
		df['zscore']= df['numerator'] / df['denominator']**(1/2)
		df['pvalue']= 2*st.norm.cdf(-abs(df['zscore']))
		df['nonmissing']= df[cols].apply(lambda x: x.count(), axis=1)
		df= df.loc[df.nonmissing >= 5, :]
		df['pos1']= np.where(pd.notna(df.pos10), df.pos10, np.where(pd.notna(df.pos11), df.pos11, df.pos12))
		df['pos2']= np.where(pd.notna(df.pos20), df.pos20, np.where(pd.notna(df.pos21), df.pos21, df.pos22))
		df['inorout']= multi.multipletests(df.pvalue, alpha= 0.05, method= 'fdr_bh')[0]
		dfHC= df.loc[df.inorout== True, :]
		df['inorout']= multi.multipletests(df.pvalue, alpha= 0.1, method= 'fdr_bh')[0]
		dfLC= df.loc[df.inorout== True,:]
		dfHC.sort_values(by=['pvalue'], inplace= True, ascending= True)
		dfLC.sort_values(by=['pvalue'], inplace= True, ascending= True)
		dfLC= dfLC.iloc[len(dfHC.index):, :]
		dfNC= df.loc[df.inorout== False, :]
		dfHC.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['segment', 'pos1', 'pos2', 'zscore', 'pvalue'])
		dfLC.to_csv(output[1], sep= '\t', header= True, index= False, columns= ['segment', 'pos1', 'pos2', 'zscore', 'pvalue'])
		dfNC.to_csv(output[2], sep= '\t', header= True, index= False, columns= ['segment', 'pos1', 'pos2', 'zscore', 'pvalue'])

rule extract_HC:
	'List of variants for extracting genotype.'
	input:
		'/mnt/work/pol/ROH/results/HC_{sample}_cox_spont'
	output:
		temp('/mnt/work/pol/ROH/results/HC_toextract_{sample}')
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		d[['chr', 'cM1', 'cM2']]= d['segment'].str.split(':',expand=True)
		d['chr']= d.chr.astype(float)
		d= d[['chr', 'pos1', 'pos2']]
		d= d.applymap(np.int64)
		d.to_csv(output[0], sep= '\t', index= False, header= False)

rule extract_vcf_samples:
        'Extract samples id included in the VCF file, for each batch.'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m12/1.vcf.gz',
                '/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m24/1.vcf.gz',
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/imputed/1.vcf.gz',
                '/mnt/archive/ROTTERDAM2/delivery-fhi/data/imputed/1.vcf.gz',
                '/mnt/archive/NORMENT1/delivery-fhi/data/imputed/feb18/1.vcf.gz',
                '/mnt/archive/NORMENT1/delivery-fhi/data/imputed/may16/1.vcf.gz'
        output:
                temp('/mnt/work/pol/ROH/{cohort}/genotypes/vcf_ids')
        run:
                if 'harvestm12' == wildcards.cohort: vcf= input[0]
                if 'harvestm24' == wildcards.cohort: vcf= input[1]
                if 'rotterdam1' == wildcards.cohort: vcf= input[2]
                if 'rotterdam2' == wildcards.cohort: vcf= input[3]
                if 'normentfeb' == wildcards.cohort: vcf= input[4]
                if 'normentmay' == wildcards.cohort: vcf= input[5]
                shell("set +o pipefail; zgrep -v '##' {vcf} | head -1 | cut -f10- | sed 's/\\t/\\n/g'  > {output[0]} ")


rule extract_samples:
        'Samples for filtering VCF files.'
        input:
                '/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv',
                '/mnt/work/pol/ROH/{cohort}/genotypes/vcf_ids'
        output:
                temp('/mnt/work/pol/ROH/{cohort}/genotypes/{sample}_ids_toextract')
        run:
                if 'harvest' in wildcards.cohort:
                        d= pd.read_csv(input[0], delim_whitespace= True, header= 0)
                        Sentrix= 'SentrixID_1'
                if 'harvest' not in wildcards.cohort:
                        d= pd.read_csv(input[0], delim_whitespace= True, header= 0)
                        Sentrix= 'SentrixID'
                if 'maternal' in wildcards.sample:
                        d= d.loc[d.Role=='Mother', :]
                if 'paternal' in wildcards.sample:
                        d= d.loc[d.Role=='Father', :]
                if 'fetal' in wildcards.sample:
                        d= d.loc[d.Role=='Child', :]
                x= [line.strip() for line in open(input[1], 'r')]
                d= d.loc[d[Sentrix].isin(x)]
                d.drop_duplicates(subset= [Sentrix], inplace= True)
                d.to_csv(output[0], header= False, columns= [Sentrix], index= False, sep= '\t')

rule extract_GT:
	'Extract genotype for HC segments.'
	input:
		'/mnt/work/pol/ROH/results/HC_toextract_{sample}',
                '/mnt/work/pol/ROH/{cohort}/genotypes/{sample}_ids_toextract',
                '/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m12/{CHR}.vcf.gz',
                '/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m24/{CHR}.vcf.gz',
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/imputed/{CHR}.vcf.gz',
                '/mnt/archive/ROTTERDAM2/delivery-fhi/data/imputed/{CHR}.vcf.gz',
                '/mnt/archive/NORMENT1/delivery-fhi/data/imputed/feb18/{CHR}.vcf.gz',
                '/mnt/archive/NORMENT1/delivery-fhi/data/imputed/may16/{CHR}.vcf.gz'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/genotypes/GT/{sample}_gt{CHR}_HC')
	run:
		if 'harvestm12' in wildcards.cohort: vcf= input[2]
		if 'harvestm24' in wildcards.cohort: vcf= input[3]
		if 'rotterdam1' in wildcards.cohort: vcf= input[4]
		if 'rotterdam2' in wildcards.cohort: vcf= input[5]
		if 'normentfeb' in wildcards.cohort: vcf= input[6]
		if 'normentmay' in wildcards.cohort: vcf= input[7]
		shell("~/soft/bcftools-1.9/bin/bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {vcf} -o {output[0]}")

rule cox_imputed:
	'Cox regression for imputed variants within HC segments.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/GT/{sample}_gt{CHR}_HC',
		'/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt',
		'/mnt/work/pol/ROH/{cohort}/genotypes/{sample}_ids_toextract'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/results/imputed/imputed_cox_spont_{sample}_temp_{CHR}')
	script:
		'scripts/cox_imputed.R'

rule cat_cox_imputed:
	'Concat results from all CHR for imputed variants.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/results/imputed/imputed_cox_spont_{{sample}}_temp_{CHR}', CHR= CHR_nms)
	output:
		'/mnt/work/pol/ROH/{cohort}/results/imputed/imputed_cox_result_{sample}'
	shell:
		'cat {input} > {output[0]}'

rule preliminary_report:
        'Generate report for harvest analysis.'
	input:
                expand('/mnt/work/pol/ROH/{{cohort}}/pheno/runs_mfr_{sample}.txt', sample= smpl_nms),
                '/mnt/work/pol/{cohort}/pheno/q1_v9.txt',
                '/mnt/work/pol/{cohort}/pheno/flag_list.txt',
                expand('/mnt/work/pol/ROH/{{cohort}}/runs/{{cohort}}_{sample}.hom', sample= smpl_nms),
                expand('/mnt/work/pol/ROH/{{cohort}}/runs/frequency/ROH_frequency_{sample}', sample= smpl_nms),
                expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{sample}.bim', pruning= pruning_nms, sample= smpl_nms),
                expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/cm_pruned{{cohort}}_{sample}.bim', pruning= pruning_nms, sample= smpl_nms),
                '/mnt/work/pol/ROH/{cohort}/runs/maternal_input_ROH_geno.txt',
                '/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',
                '/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt',
                expand('/mnt/work/pol/ROH/{{cohort}}/results/het/{sample}_excess_hom.txt', sample= smpl_nms),
                '/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt',
                '/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt',
                expand('/mnt/work/pol/ROH/results/{conf}_{sample}_cox_spont', sample= smpl_nms, conf= ['HC', 'LC', 'NC']),
		expand('/mnt/work/pol/ROH/{cohort}/results/imputed/imputed_cox_result_{sample}', cohort= cohort_nms, sample= smpl_nms),
#		expand('/mnt/work/pol/ROH/{cohort}/results/imputed/cox_spont_{sample}_{CHR}_temp', cohort= cohort_nms, sample= smpl_nms, CHR= CHR_nms),
#		expand('/mnt/work/pol/ROH/{cohort}/genotypes/GT/{sample}_gt{CHR}_HC', cohort= cohort_nms, sample= smpl_nms, CHR= CHR_nms)
        output:
                'reports/ROH_{cohort}_analysis.html'
        script:
                'scripts/report_ROH.Rmd'


rule html_meta_report:
        'Generate report for all cohorts together.'
        input:
                expand('/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt', cohort= cohort_nms, sample= smpl_nms),
                expand('/mnt/work/pol/{cohort}/pheno/q1_v9.txt', cohort= cohort_nms), 
                expand('/mnt/work/pol/ROH/{cohort}/runs/frequency/ROH_frequency_{sample}', cohort= cohort_nms, sample= smpl_nms),
                expand('/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/pruned{cohort}_{sample}.bim', cohort= cohort_nms, pruning= pruning_nms, sample= smpl_nms),
                expand('/mnt/work/pol/ROH/{cohort}/genotypes/{pruning}/cm_pruned{cohort}_{sample}.bim', cohort= cohort_nms, pruning= pruning_nms, sample= smpl_nms),
                expand('/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',cohort= cohort_nms),
                expand('/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt', cohort= cohort_nms),
                expand('/mnt/work/pol/ROH/{cohort}/results/het/{sample}_excess_hom.txt', cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt', cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt', cohort= cohort_nms),
		expand('/mnt/work/pol/{cohort}/pca/{cohort}_pca.txt', cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/{cohort}/ibd/to_phase.fam', cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_fetal.fam', cohort= cohort_nms),
		expand('/mnt/work/pol/{cohort}/pheno/flag_list.txt', cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt', cohort= cohort_nms),
		expand('/mnt/work/pol/{cohort}/relatedness/all_{cohort}.kin0', cohort= cohort_nms)
        output:
                'reports/html_meta_ROH_analysis.html'
        script:
                'scripts/html_meta_ROH.Rmd'



import pandas as pd
import numpy as np
import os
import gzip
import functools

cohort_nms= ['harvest','rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay']
smpl_nms= ['maternal','paternal', 'fetal']
batch_nms= ['m12', 'm24']
CHR_nms= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

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
		expand('/mnt/work/pol/ROH/harvest/ibd/harvest_ibd_chr{CHR}.match', CHR= CHR_nms),
		expand('/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt', cohort= cohort_nms),
		expand('reports/ROH_{cohort}_analysis.html', cohort= cohort_nms),
		expand('figures/Figure1_{cohort}.eps', cohort= cohort_nms)

## Snakemake code

rule ids_to_keep:
        'List of maternal, paternal and fetal ids acceptable by PLINK for --keep.'
        input:
                '/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt',
		'/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv'
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/maternal_ids',
                '/mnt/work/pol/ROH/{cohort}/pheno/paternal_ids',
                '/mnt/work/pol/ROH/{cohort}/pheno/fetal_ids'
	run:
		if 'harvest' in input[1]:
			d= pd.read_csv(input[0], sep= '\t')
			mat= d.loc[:,['Mother', 'Mother']]
			fet= d.loc[:, ['Child', 'Child']]
			fat= d.loc[:, ['Father', 'Father']]
	                mat.columns= ['FID', 'IID']
		        fet.columns= ['FID', 'IID']
			fat.columns= ['FID', 'IID']
		if (('rotterdam' in input[1]) | ('norment' in input[1])):
			x= pd.read_csv(input[1], delim_whitespace= True)
			x.dropna(subset= ['Role'], inplace= True)
			x.rename({'SentrixID': 'IID', 'postFID': 'FID'}, inplace= True, axis= 1)
			d= pd.read_csv(input[0], sep= '\t')
			mat= d.loc[:,['Mother', 'Mother']]
			fet= d.loc[:, ['Child', 'Child']]
			fat= d.loc[:, ['Father', 'Father']]
			mat.columns= ['Mother', 'IID']
			fet.columns= ['Child', 'IID']
			fat.columns= ['Father', 'IID']
			mat= pd.merge(mat, x, on= 'IID')
			fet= pd.merge(fet, x, on= 'IID')
			fat= pd.merge(fat, x, on= 'IID')
		mat.to_csv(output[0], header= None, columns= ['FID', 'IID'], index= False, sep= '\t')
                fet.to_csv(output[2], header= None, columns= ['FID', 'IID'], index= False, sep= '\t')
                fat.to_csv(output[1], header= None, columns= ['FID', 'IID'], index= False, sep= '\t')

rule overlaping_variants:
        'List overlapping variants between m12 and m24.'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/m12/m12-genotyped.bim',
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/m24/m24-genotyped.bim'
        output:
                temp('/mnt/work/pol/ROH/harvest/genotypes/temp/geno_to_extract')
        run:
                d12= pd.read_csv(input[0], header= None, sep= '\t')
                d24= pd.read_csv(input[1], header= None, sep= '\t')
                d12.columns= d24.columns= ['CHR','SNP','cM', 'POS', 'REF', 'ALT']
                d24.loc[d24.ALT< d24.REF, ['REF', 'ALT']]= d24.loc[d24.ALT< d24.REF, ['ALT', 'REF']]
                d12.loc[d12.ALT< d12.REF, ['REF', 'ALT']]= d12.loc[d12.ALT< d12.REF, ['ALT', 'REF']]
                d12['variant']= d12.iloc[:,0].astype(str) + ':' + d12.iloc[:,3].astype(str)+ ':' + d12.iloc[:,4] + ':' + d12.iloc[:,5]
                d24['variant']= d24.iloc[:,0].astype(str) + ':' + d24.iloc[:,3].astype(str)+ ':' + d24.iloc[:,4] + ':' + d24.iloc[:,5]
                d= pd.merge(d12, d24, on= ['variant', 'CHR', 'POS'])
                d= d.loc[:, ['CHR', 'POS', 'POS', 'variant']]
                d= d.loc[d.CHR !=23 , :]
                d.to_csv(output[0], header= False, sep= '\t', index= False)

rule overlapping_plink_geno:
        'Extract overlapping SNPs from each batch.'
        input:
                expand('/mnt/archive/HARVEST/delivery-fhi/data/genotyped/{{batch}}/{{batch}}-genotyped.{ext}', ext= ['bed','bim','fam']),
                '/mnt/work/pol/ROH/harvest/genotypes/temp/geno_to_extract'
        output:
                temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/{{batch}}_genotyped.{ext}', ext= ['bed','bim','fam', 'log']))
        params:
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/{batch}/{batch}-genotyped',
                '/mnt/work/pol/ROH/harvest/genotypes/temp/{batch}_genotyped'
        shell:
                '~/soft/plink --bfile {params[0]} --extract range {input[3]} --make-bed --out {params[1]}'

rule merge_batches_harvest_genotyped:
        'Merge batch PLINK files for ROH calling.'
        input:
                expand('/mnt/work/pol/ROH/harvest/genotypes/temp/{batch}_genotyped.{ext}', batch= batch_nms, ext= ['bed','bim','fam']),
        output:
                temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_genotyped.{ext}', ext= ['bed','bim','fam', 'log', 'nosex']))
        params:
                '/mnt/work/pol/ROH/harvest/genotypes/temp/m12_genotyped',
                '/mnt/work/pol/ROH/harvest/genotypes/temp/m24_genotyped',
                '/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_genotyped'
        shell:
                '~/soft/plink --bfile {params[0]} --bmerge {params[1]} --merge-equal-pos --make-bed --out {params[2]}'

rule copy_rott_genotyped:
        'Copy genotype PLINK files and change name.'
        input:
                expand('/mnt/archive/ROTTERDAM1/delivery-fhi/data/genotyped/genotyped.{ext}', ext= ['bed', 'bim', 'fam']),
                expand('/mnt/archive/ROTTERDAM2/delivery-fhi/data/genotyped/genotyped.{ext}', ext= ['bed', 'bim', 'fam']),
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

rule overlaping_variants_phasing:
	'List overlapping variants between m12 and m24.'
	input:
		'/mnt/archive/HARVEST/delivery-fhi/data/to_imputation/m12/m12-ready-for-imputation.bim',
		'/mnt/archive/HARVEST/delivery-fhi/data/to_imputation/m24/m24-ready-for-imputation.bim'
	output:
		temp('/mnt/work/pol/ROH/harvest/genotypes/temp/to_extract')
	run:
		d12= pd.read_csv(input[0], header= None, sep= '\t')
		d24= pd.read_csv(input[1], header= None, sep= '\t')
		d12.columns= d24.columns= ['CHR','SNP','cM', 'POS', 'REF', 'ALT']
		d24.loc[d24.ALT< d24.REF, ['REF', 'ALT']]= d24.loc[d24.ALT< d24.REF, ['ALT', 'REF']]
		d12.loc[d12.ALT< d12.REF, ['REF', 'ALT']]= d12.loc[d12.ALT< d12.REF, ['ALT', 'REF']]
		d12['variant']= d12.iloc[:,0].astype(str) + ':' + d12.iloc[:,3].astype(str)+ ':' + d12.iloc[:,4] + ':' + d12.iloc[:,5]
		d24['variant']= d24.iloc[:,0].astype(str) + ':' + d24.iloc[:,3].astype(str)+ ':' + d24.iloc[:,4] + ':' + d24.iloc[:,5]
		d= pd.merge(d12, d24, on= ['variant', 'CHR', 'POS'])
		d= d.loc[:, ['CHR', 'POS', 'POS', 'variant']]
		d= d.loc[d.CHR !=23 , :]
		d.to_csv(output[0], header= False, sep= '\t', index= False)

rule overlaping_plink_phasing:
	'Extract overlaping SNPs from each batch.'
	input:
		expand('/mnt/archive/HARVEST/delivery-fhi/data/to_imputation/{{batch}}/{{batch}}-ready-for-imputation.{ext}', ext= ['bed','bim','fam']),
		'/mnt/work/pol/ROH/harvest/genotypes/temp/to_extract'
	output:
		temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/{{batch}}_ready_for_imputation.{ext}', ext= ['bed','bim','fam', 'log']))
	params:
		'/mnt/archive/HARVEST/delivery-fhi/data/to_imputation/{batch}/{batch}-ready-for-imputation',
		'/mnt/work/pol/ROH/harvest/genotypes/temp/{batch}_ready_for_imputation'
	shell:
		'~/soft/plink --bfile {params[0]} --extract range {input[3]} --make-bed --out {params[1]}'
		
rule merge_batches_harvest_phasing:
	'Merge batch PLINK files for phasing.'
	input:
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/m12_ready_for_imputation.{ext}', ext= ['bed','bim','fam']),
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/m24_ready_for_imputation.{ext}', ext= ['bed','bim','fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing.{ext}', ext= ['bed','bim','fam', 'log', 'nosex']))
	params:
		'/mnt/work/pol/ROH/harvest/genotypes/temp/m12_ready_for_imputation',
		'/mnt/work/pol/ROH/harvest/genotypes/temp/m24_ready_for_imputation',
		'/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing'
	shell:
		'~/soft/plink --bfile {params[0]} --bmerge {params[1]} --merge-equal-pos --make-bed --out {params[2]}'

rule copy_rotterdam1_to_phasing:
        'Copy and rename PLINK files for phasing.'
        input:
                expand('/mnt/archive/ROTTERDAM1/delivery-fhi/data/to_phasing/merged/hrc-update-complete-all.{ext}', ext= ['bed', 'bim', 'fam']),
                expand('/mnt/archive/ROTTERDAM2/delivery-fhi/data/to_phasing/merged/hrc-update-complete.{ext}', ext= ['bed', 'bim', 'fam']),
		expand('/mnt/archive/NORMENT1/delivery-fhi/data/to_phasing/feb18/merged/hrc-update-complete.{ext}', ext= ['bed', 'bim', 'fam']),
		expand('/mnt/archive/NORMENT1/delivery-fhi/data/to_phasing/may16/merged/hrc-update-complete.{ext}', ext= ['bed', 'bim', 'fam'])
        output:
                temp(expand('/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_to_phasing.{ext}', ext= ['bed', 'bim', 'fam'])),
                temp(expand('/mnt/work/pol/ROH/rotterdam2/genotypes/temp/rotterdam2_to_phasing.{ext}', ext= ['bed', 'bim', 'fam'])),
        shell:
                '''
                cp {input[0]} {output[0]}; cp {input[1]} {output[1]}; cp {input[2]} {output[2]}
                cp {input[3]} {output[3]}; cp {input[4]} {output[4]}; cp {input[5]} {output[5]}
                '''

rule copy_norment_to_phasing:
        'Copy and rename PLINK files for phasing.'
        input:
                expand('/mnt/archive/NORMENT1/delivery-fhi/data/to_phasing/feb18/merged/hrc-update-complete.{ext}', ext= ['bed', 'bim', 'fam']),
                expand('/mnt/archive/NORMENT1/delivery-fhi/data/to_phasing/may16/merge/hrc-update-complete.{ext}', ext= ['bed', 'bim', 'fam'])
        output:
                temp(expand('/mnt/work/pol/ROH/normentfeb/genotypes/temp/normentfeb_to_phasing.{ext}', ext= ['bed', 'bim', 'fam'])),
                temp(expand('/mnt/work/pol/ROH/normentmay/genotypes/temp/normentmay_to_phasing.{ext}', ext= ['bed', 'bim', 'fam']))
        shell:
                '''
                cp {input[0]} {output[0]}; cp {input[1]} {output[1]}; cp {input[2]} {output[2]}
                cp {input[3]} {output[3]}; cp {input[4]} {output[4]}; cp {input[5]} {output[5]}
                '''


rule split_PLINK_chr:
	'Split PLINK binary files for phasing into one file per chromosome.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_to_phasing.{ext}', ext= ['bed','bim','fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_to_phasing_chr{{CHR}}.{ext}', ext= ['bed','bim','fam']))
	params:
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_to_phasing',
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_to_phasing_chr{CHR}'
	shell:
		'~/soft/plink --bfile {params[0]} --chr {wildcards.CHR} --mind --make-bed --out {params[1]}'

rule eagle_phasing:
	'Phasing with eagle.'
	input:
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/temp/{{cohort}}_to_phasing_chr{{CHR}}.{ext}', ext= ['bed','bim','fam'])
	output:
		temp('/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.haps.gz'),
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.sample'
	params:
		'/mnt/work/pol/ROH/{cohort}/genotypes/temp/{cohort}_to_phasing_chr{CHR}',
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}'
	shell:
		"~/soft/Eagle_v2.4.1/eagle --bfile={params[0]} --geneticMapFile={input[0]} --numThreads=10 --outPrefix={params[1]}"

rule ungzip_haps:
	'UnGzip haps output from eagle2.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.haps.gz'
	output:
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.haps'
	shell:
		'gzip -d {input[0]}'

rule haps_to_ped:
	'Format file as .ped required by GERMLINE.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.haps',
		'/mnt/work/pol/ROH/{cohort}/genotypes/haps/{cohort}_phased_chr{CHR}.sample'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_chr{CHR}.ped'),
		temp('/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_chr{CHR}.map')
	params:
		'/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_chr{CHR}'
	run:
		shell("/home/pol.sole.navais/soft/germline-1-5-3/bin/impute_to_ped {input[0]} {input[1]} {params[0]} || true")

rule add_genetic_map:
	'Adding genetic map to .map file generated.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_chr{CHR}.map',
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_complete_chr{CHR}.map')
	run:
		d= pd.read_csv(input[0], sep= ' ', header= None)
		g= pd.read_csv(input[1], sep= ' ', header= 0)
                g['SNP']= g.chr.map(str) + ':' + g.position.map(str)
                g= g[['chr', 'SNP', 'Genetic_Map(cM)', 'position']]
		d.columns= g.columns
                g= g.loc[g.chr== int(wildcards.CHR), :]
                d= d.loc[:, ['chr','position']]
                df= pd.merge(g, d, on= ['chr', 'position'], how= 'right')
                df= df.loc[df['Genetic_Map(cM)'].isna(), :]
                df['Genetic_Map(cM)']= np.interp(df.position, g['position'], g['Genetic_Map(cM)'])
                g= g.append(df)
                g['SNP']= g.chr.map(str) + ':' + g.position.map(str)
		d= pd.merge(d, g, on= ['chr', 'position'], how= 'left')
		d= d[['chr', 'SNP', 'Genetic_Map(cM)', 'position']]
                d.to_csv(output[0], header= False, index= False, sep= '\t')

rule ibd_GERMLINE:
	'Estimate shared IBD segments between subjects using GERMLINE.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_chr{CHR}.ped',
		'/mnt/work/pol/ROH/{cohort}/genotypes/ibd/{cohort}_phased_complete_chr{CHR}.map'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/ibd/{cohort}_ibd_chr{CHR}.match'),
		temp('/mnt/work/pol/ROH/{cohort}/ibd/{cohort}_ibd_chr{CHR}.log')
	params:
		'/mnt/work/pol/ROH/{cohort}/ibd/{cohort}_ibd_chr{CHR}'
	run:	
		shell("~/soft/germline-1-5-3/bin/germline -input {input[0]} {input[1]} -min_m 2 -output {params[0]} || true")

rule lightweight_ibd:
	'Remove columns from GERMLINE ibd file.'
	input:
		'/mnt/work/pol/ROH/{cohort}/ibd/{cohort}_ibd_chr{CHR}.match'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/ibd/lw_{cohort}_ibd_chr{CHR}.match')
	shell:
		"cut -d$'\t' -f1-4,7 {input[0]} > {output[0]}"

rule filter_ibd:
	'Keep only parental pairs.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/ibd/lw_{{cohort}}_ibd_chr{CHR}.match', CHR= CHR_nms),
		'/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt'
	output:
		'/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt'
	run:
		trio= [file for file in input if 'trios' in file]
		trio= pd.read_csv("".join(trio), sep= '\t', header= 0)
		df= pd.DataFrame()
		flist= [file for file in input if 'ibd' in file]
		for infile in flist:
			d= pd.read_csv(infile, header= None, delim_whitespace= True)
			d.columns= ['FID1', 'IID1', 'FID2', 'IID2', 'CHR', 'start', 'end', 'cM']
			d= d.loc[((d.IID1.isin(trio.Father.values)) & (d.IID2.isin(trio.Mother.values))) | ((d.IID1.isin(trio.Mother.values)) & (d.IID2.isin(trio.Father.values))), : ]
			d['Mother']= np.where(d.IID1.isin(trio.Mother.values), d.IID1, d.IID2)
			d['Father']= np.where(d.Mother != d.IID1, d.IID1, d.IID2)
			d= pd.merge(d, trio, on= ['Mother', 'Father'])
			d= d[['Mother','Father','Child','CHR', 'start', 'end', 'cM']]
			df= df.append(d)
		df= df.groupby(['Mother', 'Father', 'Child'])['cM'].sum().reset_index()
		df.to_csv(output[0], sep= '\t', header= True, index= False)

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
                '/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_pca.txt',
                '/mnt/work/pol/ROH/{cohort}/ibd/to_phase.fam',
                '/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_fetal.fam',
                '/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt',
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

rule combine_pca:
        'Obtain pca for all samples.'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m12-founders-pca-covariates',
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m24-founders-pca-covariates',
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m12-offspring-pca-covariates',
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m24-offspring-pca-covariates',
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/aux/pca-covar/offspring/pca/final_pca_covars.txt',
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/aux/pca-covar/founders/pca/final_pca_covars.txt'
        output:
                '/mnt/work/pol/ROH/harvest/pheno/harvest_pca.txt',
                '/mnt/work/pol/ROH/rotterdam1/pheno/rotterdam1_pca.txt'
        shell:
                '''
                cat /mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/*-pca-covariates > {output[0]}
                cat {input[4]} {input[5]} > {output[1]}
                '''

rule relatedness:
        'Calculate relatedness using KING function from PLINK2.'
        input:
                '/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_{sample}.bed',
                '/mnt/work/pol/ROH/{cohort}/pheno/{sample}_ids',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/none/pruned{{cohort}}_{{sample}}.{ext}', ext= ['bed','bim','fam']),
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}.kin0'
        params:
                '/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_{sample}',
		'/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}'
        shell:
                '~/soft/plink2 --bfile {params[0]} --keep {input[1]} --make-king-table --king-table-filter 0.03125 --out {params[1]}'

rule phenofile:
        'Merge all data necessary to create a phenotype file with ROH.'
        input:
                '/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom',
                '/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom.indiv',
                '/mnt/work/pol/{cohort}/pheno/{cohort}_mfr.csv',
                '/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv',
                '/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_pca.txt',
                '/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}.kin0',
		'/mnt/archive/HARVEST/delivery-fhi/data/genotyped/m12/m12-genotyped.fam',
		'/mnt/work/pol/ROH/{cohort}/runs/{sample}_input_ROH_geno.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{{sample}}.bim', pruning= pruning_nms)
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
        script:
                'scripts/pheno_file.py'

rule mapping_ROHs:
        'Obtain matrix (rows= position, columns = subject), with all ROHs per subject (1= homozygous part of ROH).'
        input:
                '/mnt/work/pol/ROH/{cohort}/runs/{cohort}_{sample}.hom',
		'/mnt/work/pol/ROH/{cohort}/runs/{sample}_input_ROH_geno.txt' 
        output:
                temp('/mnt/work/pol/ROH/{cohort}/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz')
        script:
                'scripts/map_ROHs.py'

rule ROH_freq:
        'Count per-position relative frequency of ROHs.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/maps/{{sample}}/maps_{{sample}}_chr{CHR}.txt.gz', CHR= CHR_nms)
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
		'/mnt/work/pol/ROH/{cohort}/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz',
		'/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
	output:
		'/mnt/work/pol/ROH/{cohort}/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}'
	script:
		'scripts/cox_ROH.R'

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
		d= functools.reduce(lambda x, y: pd.merge(x, y, on= 'IID'), dflist)
		d.to_csv(output[0], index=False, header= True, sep= '\t')

rule preliminary_report:
        'Generate report for harvest analysis.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/pheno/runs_mfr_{sample}.txt', sample= smpl_nms),
                '/mnt/work/pol/harvest/pheno/q1_pdb1724_v9.csv',
		'/mnt/work/pol/rotterdam1/pheno/q1_pdb315_v9.csv',
                expand('/mnt/work/pol/ROH/{{cohort}}/runs/{{cohort}}_{sample}.hom', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{{cohort}}/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}', CHR= CHR_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{{cohort}}/runs/frequency/ROH_frequency_{sample}', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{sample}.bim', pruning= pruning_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/cm_pruned{{cohort}}_{sample}.bim', pruning= pruning_nms, sample= smpl_nms),
		'/mnt/work/pol/ROH/{cohort}/runs/maternal_input_ROH_geno.txt',
		'/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt',
		'/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/results/het/{sample}_excess_hom.txt', sample= smpl_nms),
		'/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt',
		'/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt'
        output:
                'reports/ROH_{cohort}_analysis.html'
	script:
		'scripts/report_ROH.Rmd'

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
		if wildcards.cohort== 'harvest':
	                d12= pd.read_csv(input[0], delim_whitespace= True, header=None)
		        d24= pd.read_csv(input[1], delim_whitespace= True, header=None)
			d= pd.concat([d12, d24])
		if wildcards.cohort== 'rotterdam1':
			d= pd.read_csv(input[2], delim_whitespace= True, header= None)
		if wildcards.cohort== 'rotterdam2':
			d= pd.read_csv(input[3], delim_whitespace= True, header= None)
		if wildcards.cohort== 'normentfeb':
			d= pd.read_csv(input[4], delim_whitespace= True, header= None)
		if wildcards.cohort== 'normentmay':
			d= pd.read_csv(input[5], delim_whitespace= True, header= None)
		d.to_csv(output[0], sep= '\t', index= False, header= False)


rule figure1:
	'Figure 1 for results section.'
	input:
		'/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_pca.txt',
		'/mnt/work/pol/ROH/{cohort}/ibd/to_phase.fam',
		'/mnt/work/pol/ROH/{cohort}/genotypes/none/pruned{cohort}_fetal.fam',
		'/mnt/work/pol/ROH/{cohort}/ibd/parental_ibd.txt',
		'/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_fetal.txt'
	output:
		'figures/Figure1_{cohort}.eps'
	script:
		'scripts/figure1.R'


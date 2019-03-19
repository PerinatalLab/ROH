import pandas as pd
import numpy as np
import os
import gzip

cohort_nms= ['harvest','rotterdam1']
smpl_nms= ['maternal','paternal', 'fetal']
batch_nms= ['m12', 'm24']
CHR_nms= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

# Other arguments:
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

rule all:
	'Collect the main outputs of the workflow.'
	input:
		expand('/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt', cohort= cohort_nms, sample= smpl_nms),
#		expand('/mnt/work/pol/ROH/{cohort}/genotypes/maps/gene/{cohort}_{sample}_CHR{CHR}', cohort= cohort_nms, sample= smpl_nms, CHR= CHR_nms),
#		expand('/mnt/work/pol/ROH/{cohort}/results/maps_cox/gene/cox_spont{sample}',cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/harvest/ibd/harvest_ibd_chr{CHR}.match', CHR= CHR_nms),
		'/mnt/work/pol/ROH/arguments/arg_R2.txt',
		'/mnt/work/pol/ROH/arguments/max_R2.txt'

rule exclude_multi_allelic_rott:
	'Set range file for multi-allelic SNP detected in ROTTERDAM1.'
	output:
		temp('/mnt/work/pol/ROH/rotterdam1/multiallelic.txt')
	shell:
		'''
		echo "14 24681025        24681025        Multi" > {output}
		'''

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
		if 'rotterdam1' in input[1]:
			x= pd.read_csv(input[1], sep= ' ')
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

rule harvest_plink_split_bed:
        'Modify the bed file: remove CHR 23, 24, 25 and 26, maf <=0.05 and split file by sample. (HARVEST)'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/{batch}/{batch}-genotyped.bed',
                '/mnt/work/pol/ROH/harvest/pheno/{sample}_ids'
        output:
                temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest{{batch}}_{{sample}}.{ext}', ext= ['bed','bim','fam','prune.out','prune.in', 'log']))
        params:
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/{batch}/{batch}-genotyped',
                '/mnt/work/pol/ROH/harvest/genotypes/temp/harvest{batch}_{sample}'
        shell:
                '~/soft/plink --bfile {params[0]} --indep-pairwise 50 5 0.5 --maf 0.05 --keep {input[1]} --make-bed --not-chr 23,24,25,26 --make-founders --out {params[1]}'

rule rott_plink_split_bed:
        'Modify the bed file: remove CHR 23, 24, 25 and 26, maf <=0.05 and split file by sample. (ROTTERDAM1)'
        input:
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/genotyped/genotyped.bed',
                '/mnt/work/pol/ROH/rotterdam1/pheno/{sample}_ids',
		'/mnt/work/pol/ROH/rotterdam1/multiallelic.txt'
        output:
                temp(expand('/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_{{sample}}.{ext}', ext= ['bed','bim','fam','prune.out','prune.in', 'log']))
        params:
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/genotyped/genotyped',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_{sample}'
        shell:
                '~/soft/plink --bfile {params[0]} --exclude range {input[2]} --indep-pairwise 50 5 0.5 --maf 0.05 --keep {input[1]} --make-bed --not-chr 23,24,25,26 --make-founders --out {params[1]}'

rule harvest_plink_bfile_prune:
        'Exclude genetic variants in prune.out files (obtained with rule plink_split_bed). (HARVEST)'
        input:
                expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest{{batch}}_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'prune.out'])
        output:
                expand('/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{{batch}}_{{sample}}.{ext}', ext= ['bed','bim','fam'])
        params:
                '/mnt/work/pol/ROH/harvest/genotypes/temp/harvest{batch}_{sample}',
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}',
                '/mnt/work/pol/ROH/harvest/genotypes/temp/harvest{batch}_{sample}.prune.out'
        shell:
                '~/soft/plink --bfile {params[0]} --exclude {params[2]} --make-bed --out {params[1]}'

rule rotterdam1_plink_bfile_prune:
        'Exclude genetic variants in prune.out files (obtained with rule plink_split_bed). (ROTTERDAM1)'
        input:
                expand('/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam', 'prune.out'])
        output:
                expand('/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{{sample}}.{ext}', ext= ['bed','bim','fam'])
        params:
                '/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_{sample}',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/temp/rotterdam1_{sample}.prune.out'
        shell:
                '~/soft/plink --bfile {params[0]} --exclude {params[2]} --make-bed --out {params[1]}'

rule overlaping_variants:
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

rule overlaping_plink:
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

rule split_PLINK_chr:
	'Split PLINK binary files for phasing into one file per chromosome.'
	input:
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing.{ext}', ext= ['bed','bim','fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing_chr{{CHR}}.{ext}', ext= ['bed','bim','fam']))
	params:
		'/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing',
		'/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing_chr{CHR}'
	shell:
		'~/soft/plink --bfile {params[0]} --chr {wildcards.CHR} --mind --make-bed --out {params[1]}'

rule eagle_phasing:
	'Phasing with eagle.'
	input:
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt',
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing_chr{{CHR}}.{ext}', ext= ['bed','bim','fam'])
	output:
		temp('/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.haps.gz'),
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.sample'
	params:
		'/mnt/work/pol/ROH/harvest/genotypes/temp/harvest_to_phasing_chr{CHR}',
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}'
	shell:
		"~/soft/Eagle_v2.4.1/eagle --bfile={params[0]} --geneticMapFile={input[0]} --numThreads=10 --outPrefix={params[1]}"

rule ungzip_haps:
	'UnGzip haps output from eagle2.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.haps.gz'
	output:
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.haps'
	shell:
		'gzip -d {input[0]}'

rule haps_to_ped:
	'Format file as .ped required by GERMLINE.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.haps',
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.sample'
	output:
		temp('/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_chr{CHR}.ped'),
		temp('/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_chr{CHR}.map')
	params:
		'/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_chr{CHR}'
	run:
		shell("/home/pol.sole.navais/soft/germline-1-5-3/bin/impute_to_ped {input[0]} {input[1]} {params[0]} || true")

rule add_genetic_map:
	'Adding genetic map to .map file generated.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_chr{CHR}.map',
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt'
	output:
		temp('/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_complete_chr{CHR}.map')
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
		'/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_chr{CHR}.ped',
		'/mnt/work/pol/ROH/harvest/genotypes/ibd/harvest_phased_complete_chr{CHR}.map'
	output:
		'/mnt/work/pol/ROH/harvest/ibd/harvest_ibd_chr{CHR}.match'
	params:
		'/mnt/work/pol/ROH/harvest/ibd/harvest_ibd_chr{CHR}'
	run:	
		shell("~/soft/germline-1-5-3/bin/germline -input {input[0]} {input[1]} -output {params[0]} || true")

rule lightweight_ibd:
	'Remove columns from GERMLINE ibd file.'
	input:
		'/mnt/work/pol/ROH/harvest/ibd/harvest_ibd_chr{CHR}.match'
	output:
		temp('/mnt/work/pol/ROH/harvest/ibd/lw_harvest_ibd_chr{CHR}.match')
	shell:
		"cut -d$'\t' -f1-4,7 {input[0]} > {output[0]}"

rule filter_ibd:
	'Keep only parental pairs.'
	input:
		expand('/mnt/work/pol/ROH/harvest/ibd/lw_harvest_ibd_chr{CHR}.match', CHR= CHR_nms),
		'/mnt/work/pol/ROH/harvest/pheno/harvest_trios.txt'
	output:
		'/mnt/work/pol/ROH/harvest/ibd/parental_ibd.txt'
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
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.bim',
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt'
	output:
		'/mnt/work/pol/ROH/harvest/genotypes/cm_prunedharvest_{batch}_{sample}.bim'
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
		df['pos']= df['X']
		df= df[['chr', 'SNP', 'X', 'pos', 'A1', 'A2']]
		df.to_csv(output[0], sep= '\t', header= False, index= False)

rule run_ROH_multi_arg:
	'Estimate ROH using multiple arguments.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_fetal.bed',
		'/mnt/work/pol/ROH/harvest/genotypes/cm_prunedharvest_{batch}_fetal.bim',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_fetal.fam'
	output:
		temp(expand('/mnt/work/pol/ROH/harvest/multi/fetal{{dens}}_{{SNP}}_{{length}}_{{het}}_{{GAP}}_{{batch}}.{ext}', ext= ['log', 'hom', 'hom.indiv', 'nosex', 'hom.summary']))
	params:
		'/mnt/work/pol/ROH/harvest/multi/fetal{dens}_{SNP}_{length}_{het}_{GAP}_{batch}'
	shell:
		'''
		/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --homozyg-window-snp {wildcards.SNP}*0.05 --homozyg-snp {wildcards.SNP} --homozyg-kb {wildcards.length} --homozyg-gap {wildcards.GAP} --homozyg-window-missing {wildcards.SNP}*0.0005 --homozyg-window-het {wildcards.het} --homozyg-density {wildcards.dens} --out {params[0]}
		'''

rule run_ROH_multi_arg_bp:
        'Estimate ROH using multiple arguments.'
        input:
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_fetal.bed',
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_fetal.bim',
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_fetal.fam'
        output:
                temp(expand('/mnt/work/pol/ROH/harvest/multi/bpfetal{{densbp}}_{{SNPbp}}_{{lengthbp}}_{{hetbp}}_{{GAPbp}}_{{batch}}.{ext}', ext= ['log', 'hom', 'hom.indiv', 'nosex', 'hom.summary']))
        params:
                '/mnt/work/pol/ROH/harvest/multi/bpfetal{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}_{batch}'
        shell:
                '''
                /home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --homozyg-window-snp {wildcards.SNPbp}*0.05 --homozyg-snp {wildcards.SNPbp} --homozyg-kb {wildcards.lengthbp} --homozyg-gap {wildcards.GAPbp} --homozyg-window-missing {wildcards.SNPbp}*0.0005 --homozyg-window-het {wildcards.hetbp} --homozyg-density {wildcards.densbp} --out {params[0]}
                '''

rule merge_batch_multi_ROH:
	'Merge batches from multiple argument ROH.'
	input:
		'/mnt/work/pol/ROH/harvest/multi/fetal{dens}_{SNP}_{length}_{het}_{GAP}_m12.hom.indiv',
		'/mnt/work/pol/ROH/harvest/multi/fetal{dens}_{SNP}_{length}_{het}_{GAP}_m24.hom.indiv'
	output:
		'/mnt/work/pol/ROH/harvest/multi/fetal_allbatch_{dens}_{SNP}_{length}_{het}_{GAP}.hom.indiv'
	run:
		d12= pd.read_csv(input[0], delim_whitespace= True, header=0)
		d24= pd.read_csv(input[1], delim_whitespace= True, header=0)
		d= pd.concat([d12, d24])
		d.to_csv(output[0], index= False, header= True, sep= '\t')


rule merge_batch_multi_ROH_bp:
        'Merge batches from multiple argument ROH.'
        input:
                '/mnt/work/pol/ROH/harvest/multi/bpfetal{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}_m12.hom.indiv',
                '/mnt/work/pol/ROH/harvest/multi/bpfetal{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}_m24.hom.indiv'
        output:
                '/mnt/work/pol/ROH/harvest/multi/bp_fetal_allbatch_{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}.hom.indiv'
        run:
                d12= pd.read_csv(input[0], delim_whitespace= True, header=0)
                d24= pd.read_csv(input[1], delim_whitespace= True, header=0)
                d= pd.concat([d12, d24])
                d.to_csv(output[0], index= False, header= True, sep= '\t')


rule determine_arguments_ROH:
	'Determine ROH estimation arguments that maximise ROH - parental IBD correlation.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm12_fetal.fam',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm24_fetal.fam',
		expand('/mnt/work/pol/ROH/harvest/multi/fetal_allbatch_{dens}_{SNP}_{length}_{het}_{GAP}.hom.indiv', dens= dens_nms, SNP= SNP_nms, length= length_nms, het= het_nms, GAP= GAP_nms),
		'/mnt/work/pol/ROH/harvest/ibd/parental_ibd.txt',
		expand('/mnt/work/pol/ROH/harvest/multi/bp_fetal_allbatch_{densbp}_{SNPbp}_{lengthbp}_{hetbp}_{GAPbp}.hom.indiv', densbp= dens_bp, SNPbp= SNP_bp, lengthbp= length_bp, hetbp= het_bp, GAPbp= GAP_bp)
	output:
		'/mnt/work/pol/ROH/arguments/arg_R2.txt',
		'/mnt/work/pol/ROH/arguments/max_R2.txt'
	run:
		ibd= [file for file in input if 'ibd' in file]
		ibd= pd.read_csv("".join(ibd), sep= '\t', header=0)
		d12= pd.read_csv(input[0], header= None, sep= ' ')
		d24= pd.read_csv(input[1], header= None, sep= ' ')
		fam= d12.append(d24)
		fam.columns= ['FID', 'IID', 'x1','x2', 'x3','x4']
		fam= fam.loc[:, ['IID', 'x1']]
		infile= list()
		R2= list()
		for f in input:
			if '.hom.indiv' in f:
				infile.append(f)
				d= pd.read_csv(f, sep= '\t', header= 0)
				d= pd.merge(d, fam, how= 'outer')
				d['KB']= np.where(d.KB.isna(), 0, d.KB)
				d= pd.merge(ibd, d, left_on= ['Child'], right_on= ['IID'])
				R2.append((d[['KB', 'cM']].corr(method= 'pearson').iloc[0,1])**2)
		d= pd.DataFrame(np.column_stack([infile, R2]), columns=['file', 'R2'])
		d.to_csv(output[0], sep= '\t', header= True, index= False)
		d= d.sort_values(['R2'], ascending= False)
		winner= d.iloc[0,0]
		winner= winner.split('_')
		winner= [i.replace('.hom.indiv', '') for i in winner]
		winner= [float(s) for s in winner if isfloat(s)]	
		with open(output[1], 'a') as out:
			out.writelines('%s\n' % value for value in winner)

rule estimate_ROH:
        '''
        Obtain ROH estimates using PLINK 1.9.
        Configuration according to file "/mnt/work/pol/ROH/arguments/max_R2.txt"
        '''
        input:
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.bed',
		'/mnt/work/pol/ROH/harvest/genotypes/cm_prunedharvest_{batch}_{sample}.bim',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.fam',
		'/mnt/work/pol/ROH/arguments/max_R2.txt',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.bim'
        output:
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom.indiv',
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom'
	params:
		'/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}'
	run:
		parlist= [line.strip() for line in open(input[3], 'r')]
		parlist= [float(x) for x in parlist]
		if parlist[0] < 100:
			shell("/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --homozyg-window-snp {parlist[1]}*0.05 --homozyg-snp {parlist[1]} --homozyg-kb {parlist[2]} --homozyg-gap {parlist[4]} --homozyg-window-missing {parlist[1]}*0.0005 --homozyg-window-het {parlist[3]} --homozyg-density {parlist[0]} --out {params}")
		if parlist[0] > 100:
			shell("/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[4]} --fam {input[2]} --homozyg-window-snp {parlist[1]}*0.05 --homozyg-snp {parlist[1]} --homozyg-kb {parlist[2]} --homozyg-gap {parlist[4]} --homozyg-window-missing {parlist[1]}*0.0005 --homozyg-window-het {parlist[3]} --homozyg-density {parlist[0]} --out {params}")

rule rott_estimate_ROH:
        '''
        Obtain ROH estimates using PLINK 1.9.
	Configuration according to file "/mnt/work/pol/ROH/arguments/max_R2.txt"
        '''
        input:
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.bed',
		'/mnt/work/pol/ROH/rotterdam1/genotypes/cm_prunedrotterdam1_{sample}.bim',
		'/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.fam',
		'/mnt/work/pol/ROH/arguments/max_R2.txt',
		'/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.bim'
        output:
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom.indiv',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom'
        params:
		'/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}'
	run:
		parlist= [line.strip() for line in open(input[3], 'r')]
		parlist= [float(x) for x in parlist]
                if parlist[0] < 100:
                        shell("/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --homozyg-window-snp parlist[1]*0.05 --homozyg-snp parlist[1] --homozyg-kb 0.0000001 --homozyg-gap parlist[4] --homozyg-window-missing parlist[1]*0.0005 --homozyg-window-het parlist[3] --homozyg-density parlist[0] --out {params[0]}")
		if parlist[0] > 100:
			shell("/home/pol.sole.navais/soft/plink --bed {input[0]} --bim {input[4]} --fam {input[2]} --homozyg-window-snp parlist[1]*0.05 --homozyg-snp parlist[1] --homozyg-kb 0.0000001 --homozyg-gap parlist[4] --homozyg-window-missing parlist[1]*0.0005 --homozyg-window-het parlist[3] --homozyg-density parlist[0] --out {params[0]}")

rule merge_m12_m24:
	'Merge PLINK files from the two HARVEST batches for relatedness calculation.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm12_{sample}.bed',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm24_{sample}.bed',
		expand('/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{{sample}}.{ext}', batch= batch_nms, ext= ['bed', 'bim','fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_{{sample}}.{ext}', ext= ['bed','bim','fam','log']))
	params:
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm12_{sample}',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm24_{sample}',
		'/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_{sample}'
	shell:
		'~/soft/plink --bfile {params[0]} --bmerge {params[1]} --merge-equal-pos --out {params[2]}'

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
                '/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_{sample}.bed',
                '/mnt/work/pol/ROH/{cohort}/pheno/{sample}_ids',
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_{{sample}}.{ext}', ext= ['bed','bim','fam']),
		expand('/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{{sample}}.{ext}', ext= ['bed', 'bim', 'fam'])
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}.kin0'
        params:
                '/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_{sample}',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}',
                '/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}'
        shell:
                '''
		if [ {wildcards.cohort} == 'harvest' ]
                then
                        ~/soft/plink2 --bfile {params[0]} --keep {input[1]} --make-king-table --king-table-filter 0.03125 --out {params[2]}
                elif [ {wildcards.cohort} == 'rotterdam1' ]
                then
                        ~/soft/plink2 --bfile {params[1]} --keep {input[1]} --make-king-table --king-table-filter 0.03125 --out {params[2]}
                fi
                '''

rule phenofile:
        'Merge all data necessary to create a phenotype file with ROH.'
        input:
                '/mnt/work/pol/ROH/harvest/runs/harvestm12_{sample}.hom',
                '/mnt/work/pol/ROH/harvest/runs/harvestm12_{sample}.hom.indiv',
                '/mnt/work/pol/ROH/harvest/runs/harvestm24_{sample}.hom',
                '/mnt/work/pol/ROH/harvest/runs/harvestm24_{sample}.hom.indiv',
                '/mnt/work/pol/harvest/pheno/harvest_mfr.csv',
                '/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv',
                '/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_pca.txt',
                '/mnt/work/pol/ROH/harvest/genotypes/cm_prunedharvestm12_{sample}.bim',
                '/mnt/work/pol/ROH/harvest/genotypes/cm_prunedharvestm24_{sample}.bim',
                '/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}.kin0',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom.indiv',
                '/mnt/work/pol/rotterdam1/pheno/rotterdam1_MFR.csv',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/cm_prunedrotterdam1_{sample}.bim',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm12_{sample}.fam'
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
        script:
                'scripts/pheno_file.py'

rule mapping_ROHs:
        'Obtain matrix (rows= position, columns = subject), with all ROHs per subject (1= homozygous part of ROH).'
        input:
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom',
                '/mnt/work/pol/ROH/harvest/genotypes/cm_prunedharvest{batch}_{sample}.bim',
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.fam'
        output:
                temp('/mnt/work/pol/ROH/harvest/genotypes/maps/{sample}/maps_{batch}_{sample}_chr{CHR}.txt')
        params:
                '/mnt/work/pol/ROH/harvest/genotypes/maps/{sample}/maps_{batch}_{sample}_chr'
        script:
                'scripts/map_ROHs.py'

rule rott_mapping_ROHs:
        'Obtain matrix (rows= position, columns = subject), with all ROHs per subject (1= homozygous part of ROH).'
        input:
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/cm_prunedrotterdam1_{sample}.bim',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.fam'
        output:
                '/mnt/work/pol/ROH/rotterdam1/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt'
        params:
                '/mnt/work/pol/ROH/rotterdam1/genotypes/maps/{sample}/maps_{sample}_chr'
        script:
                'scripts/map_ROHs.py'

rule merge_maps:
        'Merge maps from the two batches, one file per chromosome and sample.'
        input:
                '/mnt/work/pol/ROH/harvest/genotypes/maps/{sample}/maps_m12_{sample}_chr{CHR}.txt',
                '/mnt/work/pol/ROH/harvest/genotypes/maps/{sample}/maps_m24_{sample}_chr{CHR}.txt'
        output:
                '/mnt/work/pol/ROH/harvest/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz'
        run:
                d12= pd.read_csv(input[0], sep= '\t', header= 0)
                d24= pd.read_csv(input[1], sep= '\t', header= 0)
                d= pd.merge(d12, d24, how= 'outer', on =['CHR', 'BP'])
                d.to_csv(output[0], sep= '\t', compression= 'gzip', index= False)

rule gzip_rott_maps:
	'Gzip rotterdam ROH maps.'
	input:
		'/mnt/work/pol/ROH/rotterdam1/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt'
	output:
		'/mnt/work/pol/ROH/rotterdam1/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz'
	shell:
		'gzip {input}'

rule ROH_freq:
        'Count per-position relative frequency of ROHs.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/maps/{{sample}}/maps_{{sample}}_chr{CHR}.txt.gz', CHR= CHR_nms)
        output:
                '/mnt/work/pol/ROH/{cohort}/runs/frequency/ROH_frequency_{sample}'
        run:
                for i in input:
                        for chunk in pd.read_csv(gzip.open(i), sep ='\t', index_col= 0, chunksize= 500):
                                chunk.fillna(0, inplace= True)
                                x= chunk.iloc[:,2:].mean(axis= 1)
                                x= pd.concat([chunk.iloc[:,0:2], x], axis= 1, ignore_index= True, sort= False)
                                x.to_csv(output[0], mode= 'a', sep= '\t', header= False, index= False)

rule cox_ph_analysis:
	''
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz',
		'/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
	output:
		'/mnt/work/pol/ROH/{cohort}/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}'
	script:
		'scripts/cox_ROH.R'

rule clean_UCSC_gene:
	'Clean UCSC gene list.'
	input:
		'raw_data/UCSC_hg19_gene_list'
	output:
		'raw_data/UCSC_hg19_gene_list_clean'
	run:
		df= pd.read_csv(input[0], sep= '\t')
		df.columns= ['gene', 'geneID', 'chr', 'start', 'end', 'cds', 'cde', 'ID']
		df= df.loc[df.cds != df.cde, :]
		df.drop(['geneID', 'ID', 'cds', 'cde'], axis= 1, inplace= True)
		df['chr']= df.chr.str.replace('chr','')
		df= df[df.chr.apply(lambda x: x.isnumeric())]
		df['length']= df.end - df.start
		df.sort_values(by= ['length'], ascending= False, inplace= True)
		df.drop_duplicates(['gene'], keep= 'first', inplace= True)
		df.drop(['length'], axis=1, inplace= True)
		df.to_csv(output[0], sep= '\t', index= False, header= True)

rule rott_gene_based_map:
	''
	input:
		'raw_data/UCSC_hg19_gene_list_clean',
		'/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom',
		'/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.fam'
	output:
		'/mnt/work/pol/ROH/rotterdam1/genotypes/maps/gene/rotterdam1_{sample}_CHR{CHR}'
	script:
		'scripts/gene_maps_ROH.py'

rule gene_based_map:
	''
	input:
		'raw_data/UCSC_hg19_gene_list_clean',
		'/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.fam'
	output:
		temp('/mnt/work/pol/ROH/harvest/genotypes/maps/gene/harvest{batch}_{sample}_CHR{CHR}')
	script:
		'scripts/gene_maps_ROH.py'

rule merge_gene_maps:
	'Merge genetic maps from the two batches, one file per chromosome and sample.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/maps/gene/harvestm12_{sample}_CHR{CHR}',
		'/mnt/work/pol/ROH/harvest/genotypes/maps/gene/harvestm24_{sample}_CHR{CHR}'
	output:
		'/mnt/work/pol/ROH/harvest/genotypes/maps/gene/harvest_{sample}_CHR{CHR}'
	run:
		d12= pd.read_csv(input[0], sep= '\t', header= 0)
		d24= pd.read_csv(input[1], sep= '\t', header= 0)
		d= pd.merge(d12, d24, how= 'outer', on =['chr', 'gene'])
		d.to_csv(output[0], sep= '\t', index= False) 

rule gene_ROH_cox:
	'Survivanl analysis for gene FROH on spontaneous delivery risk.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/maps/gene/{cohort}_{sample}_CHR{CHR}',
		'/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
	output:
		temp('/mnt/work/pol/ROH/{cohort}/results/maps_cox/gene/{sample}/cox_spont_{sample}_CHR{CHR}')
	script:
		'scripts/cox_gene_ROH.R'


rule trios_list:
	'Obtain a list of family trio IDs.'
	input:
		'/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv'
	output:
		'/mnt/work/pol/ROH/{cohort}/pheno/{cohort}_trios.txt'
	run:
		if 'harvest' in input[0]:
			d= pd.read_csv(input[0], sep= ';')
			d.dropna(subset= ['Role'], inplace= True)
			d= d.pivot(index= 'PREG_ID_1724', columns= 'Role', values= 'SentrixID_1')
		
		if 'rotterdam1' in input[0]:
			d= pd.read_csv(input[0], sep= ' ')
			d.dropna(subset= ['Role'], inplace= True)
			d= d.pivot(index= 'PREG_ID_315', columns= 'Role', values= 'SentrixID')
		d.dropna(inplace= True, axis= 0)
		d.reset_index(inplace= True)
		d.to_csv(output[0], header=True, sep= '\t', index= False)

rule cat_gene_based_results:
	'Concatenate all result files for each sample and cohort.'
	input:
		expand('/mnt/work/pol/ROH/{{cohort}}/results/maps_cox/gene/{{sample}}/cox_spont_{{sample}}_CHR{CHR}', CHR= CHR_nms)
	output:
		'/mnt/work/pol/ROH/{cohort}/results/maps_cox/gene/cox_spont{sample}'
	shell:
		'cat {input} > {output}'

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
		'/mnt/work/pol/ROH/harvest/genotypes/haps/harvest_phased_chr{CHR}.haps'
	output:
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37_{CHR}.map'
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

rule generate_report:
        'Generate report for harvest analysis.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/pheno/runs_mfr_{sample}.txt', sample= smpl_nms),
                '/mnt/work/pol/{cohort}/pheno/mod_{cohort}_q1_v9.csv',
                expand('/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom', sample= smpl_nms, batch= batch_nms),
                expand('/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom', sample= smpl_nms),
#                expand('/mnt/work/pol/ROH/rotterdam1/genotypes/maps/gene/rotterdam1_{sample}_CHR{CHR}', sample= smpl_nms, CHR= CHR_nms),
#                expand('/mnt/work/pol/ROH/harvest/genotypes/maps/gene/harvest{batch}_{sample}_CHR{CHR}', sample= smpl_nms, batch= batch_nms, CHR= CHR_nms),
		expand('/mnt/work/pol/ROH/{{cohort}}/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}', CHR= CHR_nms, sample= smpl_nms)
        output:
                '/home/pol.sole.navais/ROH/reports/ROH_{cohort}_analysis.html'
        shell:
                """
                echo 'rmarkdown::render(input="scripts/report_ROH.Rmd", output_file="{output}")' | R --vanilla
                """



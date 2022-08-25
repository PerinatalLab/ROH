import pandas as pd
import numpy as np
import os
import gzip
from functools import reduce

cohort_nms= ['harvestm12', 'harvestm24','rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay']
smpl_nms= ['maternal','paternal', 'fetal']
batch_nms= ['m12', 'm24']
CHR_nms= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22]
rep_nms= ['normentjan', 'normentjun']
Rclass_nms= ['classA', 'classB', 'classC']

sample_rep_nms= ['maternal']
# Other arguments:

pruning_nms= ['none', 'soft', 'moderate']

dens_nms= [5]
SNP_nms= [15, 25, 50, 75, 100, 150, 200, 300, 400]
length_nms= [0.0000001]
het_nms= [0, 1]
GAP_nms= [5]

dens_bp= [5000]
SNP_bp= [15, 25, 50, 75, 100, 150, 200, 300, 400]
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
		expand('/mnt/work/pol/ROH/{cohort}/pheno/{sample}_ids.txt', cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/pheno/runs_mfr_{sample}.txt', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/arguments/arg_R2_{cohort}.txt', cohort= cohort_nms),
		expand('/mnt/work/pol/ROH/arguments/max_R2_{cohort}.txt', cohort= cohort_nms),
#		expand('/mnt/work/pol/ROH/pheno/excess_{sample}.txt', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/results/surv_spont_{sample}', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/results/imputed/surv_imputed_{sample}.txt', sample= smpl_nms),
#		'/mnt/work/pol/ROH/reports/meta_beamer_ROH.pdf',
		expand('/mnt/work/pol/ROH/annotation/independent_OMIM_HC_{sample}.txt', sample= smpl_nms),
#		'/mnt/work/pol/ROH/reports/Figures.pdf',
		'/mnt/work/pol/ROH/figures/parent_offspring_assoc_optim.eps',
		'/mnt/work/pol/ROH/figures/SNP_R2_optim.eps',
		expand('/mnt/work/pol/ROH/figures/zscore_mht_{sample}.eps', sample= smpl_nms),
		'/mnt/work/pol/ROH/figures/ROH_frequency.eps',
		'/mnt/work/pol/ROH/figures/tmrca.eps',
#	expand('/mnt/work/pol/ROH/figures/individual_segments_{sample}.eps', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/figures/segments_pvalue_{sample}.eps', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/figures/survival_curves_{sample}.eps', sample= smpl_nms),
		'/mnt/work/pol/ROH/tables/AFT_FROH.txt',
		'/mnt/work/pol/ROH/tables/descr_cohorts.txt',
		'/mnt/work/pol/ROH/tables/autoz_all.txt',
		expand('/mnt/work/pol/ROH/tables/HC_indep_annotated_{sample}.txt',sample= smpl_nms),
		expand('/mnt/work/pol/ROH/tables/HC_annotated_{sample}.txt', sample= smpl_nms),
		'/mnt/work/pol/ROH/tables/optim_param.txt',
		expand('/mnt/work/pol/ROH/figures/S{n_fig}_Figure.pdf', n_fig= [1, 2, 3, 4, 5, 6, 7, 8]),
		expand('/mnt/work/pol/ROH/tables/S{n_fig}_Table.pdf',n_fig= [1, 2, 3, 4, 5, 6, 7, 8]),
		expand('/mnt/work/pol/ROH/results/burden_survival_{sample}.txt', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/results/{sample}/gene_burden_eff_ROH.txt', sample= smpl_nms), 
		expand('/mnt/work/pol/ROH/results/{sample}/eff_ROH.txt', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/results/{sample}/loglik_{sample}.txt', sample= smpl_nms),
#		expand('/mnt/work/pol/ROH/replication/pheno/runs_mfr_{sample_rep}.txt', sample_rep= sample_rep_nms),
		expand('/mnt/work/pol/ROH/replication/results/imputed_surv_spont_maternal', sample_rep= sample_rep_nms),
		expand('/mnt/work/pol/ROH/replication/results/burden_survival_{sample_rep}.txt', sample_rep= sample_rep_nms),
		expand('/mnt/work/pol/ROH/replication/results/autoz_surv_spont_{sample_rep}', sample_rep= sample_rep_nms),
		expand('/mnt/work/pol/ROH/genotypes/lof/geno/top_missense_{sample}.txt', sample= smpl_nms),
#		expand('/mnt/work/pol/ROH/fixed_params/results/Joshi_burden_survival_{sample}.txt', sample=smpl_nms),
#		expand('/mnt/work/pol/ROH/fixed_params/results/{sample}/Joshi_eff_ROH.txt', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/fixed_params/overlap_Joshi_params_{sample}_ROH.txt', sample= smpl_nms),
		'/mnt/work/pol/ROH/figures/qqplot_segments.eps',
		expand('/mnt/work/pol/ROH/figures/gene_burden_mht_{sample}.eps', sample= smpl_nms),
		'/mnt/work/pol/ROH/figures/qqplot_gene_burden.eps',
		'/mnt/work/pol/ROH/figures/maternal_fetal_gene_burden.eps',
		'/mnt/work/pol/ROH/figures/ROH_size.eps',
		'/mnt/work/pol/ROH/figures/density_overlap.eps',
		expand('/mnt/work/pol/ROH/figures/gene_burden_pvalue_{sample}.eps', sample= smpl_nms),
		'/mnt/work/pol/ROH/figures/locus_2_gene_burden_pvalue_maternal.eps',
		'/mnt/work/pol/ROH/tables/S_Tables.pdf',
		'/mnt/work/pol/ROH/figures/S_Figures.pdf'

include: 'scripts/survival/Snakefile'
include: 'scripts/figures/Snakefile'
include: 'scripts/metaanalysis/Snakefile'
include: 'scripts/segments_snv_maps/Snakefile'
include: 'scripts/annotation/Snakefile'
include: 'scripts/reports/Snakefile'
include: 'scripts/imputed/Snakefile'
include: 'scripts/phasing/Snakefile'
include: 'scripts/ROH_calling/Snakefile'
include: 'scripts/replication/Snakefile'
include: 'scripts/chrX/Snakefile'
include: 'scripts/fixed/Snakefile'

## Snakemake code

rule ids_to_keep:
	'List maternal, paternal and fetal ids acceptable by PLINK for --keep.'
	input:
		'/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv'
	output:
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
			x.to_csv(output[3], header= True, sep= '\t', index= False)

rule ids:
	'List maternal, paternal and fetal ids acceptable by PLINK for --keep.'
	input:
		'/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv'
	output:
		'/mnt/work/pol/ROH/{cohort}/pheno/{sample}_ids.txt'
	run:
		if 'harvest' in wildcards.cohort:
			d= pd.read_csv(input[0], sep= '\t', header= 0)
			d.dropna(subset= ['Role'], inplace= True)
			x= d.pivot(index='PREG_ID_1724', columns='Role', values= 'SentrixID_1').reset_index()
			x.dropna(inplace= True)
			x_c= x.dropna(subset= ['Child'])
			x_m= x.dropna(subset= ['Mother'])
			x_f= x.dropna(subset= ['Father'])
			if wildcards.sample== 'maternal':
				x_m.to_csv(output[0], header= None, columns= ['Mother', 'Mother'], index= False, sep= '\t')
			if wildcards.sample== 'fetal':
				x_c.to_csv(output[0], header= None, columns= ['Child', 'Child'], index= False, sep= '\t')
			if wildcards.sample== 'paternal':
				x_f.to_csv(output[0], header= None, columns= ['Father', 'Father'], index= False, sep= '\t')
		if 'harvest' not in wildcards.cohort:
			d= pd.read_csv(input[0], delim_whitespace= True, header= 0)
			d.dropna(subset= ['Role'], inplace= True)
			x= d.pivot(index= 'PREG_ID_315', columns= 'Role', values= ['FID', 'SentrixID']).reset_index()
			FID= np.where(~x.FID.Child.isnull(), x.FID.Child, np.where(~x.FID.Mother.isnull(), x.FID.Mother, x.FID.Father))
			x= pd.concat([x.PREG_ID_315, x.SentrixID, pd.DataFrame({'FID': FID})], axis= 1)
			x['PREG_ID_315']= x['PREG_ID_315'].astype(int)
			x= x[pd.to_numeric(x['FID'], errors='coerce').notnull()]
			x['FID']= x['FID'].astype(int)
			x.dropna(subset= ['FID'], inplace= True)
			x_c= x.dropna(subset= ['Child'])
			x_m= x.dropna(subset= ['Mother'])
			x_f= x.dropna(subset= ['Father'])
			if wildcards.sample== 'maternal':
				x_m.to_csv(output[0], header= None, columns= ['FID', 'Mother'], index= False, sep= '\t')
			if wildcards.sample== 'fetal':
				x_c.to_csv(output[0], header= None, columns= ['FID', 'Child'], index= False, sep= '\t')
			if wildcards.sample== 'paternal':
				x_f.to_csv(output[0], header= None, columns= ['FID', 'Father'], index= False, sep= '\t')


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
		'/mnt/work/pol/{cohort}/pca/pca_exclude.txt',
		'/mnt/work/pol/ROH/{cohort}/runs/sUPD_{cohort}_{sample}.txt',
		expand('/mnt/work/pol/ROH/{{cohort}}/genotypes/{pruning}/pruned{{cohort}}_{{sample}}.bim', pruning= pruning_nms)
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
        script:
                'scripts/pheno_file.py'

rule concat_phenos_PCA:
        'Concat pheno files, and add PCA.'
        input:
                '/mnt/archive/MOBAGENETICS/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc',
                '/mnt/archive/MOBAGENETICS/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0',
                expand('/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{{sample}}.txt', cohort= cohort_nms)
        output:
                '/mnt/work/pol/ROH/pheno/runs_mfr_{sample}.txt'
        run:
                def selectUnrelated(df, x):
                        kin= pd.read_csv(input[1], header= 0, sep= '\t')
                        kin= kin.loc[kin.Kinship > 0.0884, :]
                        kin= kin.loc[kin.ID1.isin(x.values)]
                        kin= kin.loc[kin.ID2.isin(x.values)]
                        kin= kin.loc[:, ['ID1','ID2','Kinship']]
                        kin_temp= kin.copy()
                        kin_temp.columns= ['ID2', 'ID1', 'Kinship']
                        kin_temp= kin_temp.append(kin)
                        kin_temp['n']= kin_temp.groupby('ID1')['ID1'].transform('count')
                        kin_temp['nn']= kin_temp.groupby('ID2')['ID2'].transform('count')
                        kin_temp.sort_values(by=['n', 'nn'], inplace= True)
                        to_keep= list()
                        for i in range(0, len(kin_temp.index)):
                                if kin_temp.iloc[i, 0] in kin_temp.iloc[0:i, 1].values:
                                        kin_temp.iloc[i, 1]= "X"
                                else:
                                        to_keep.append(kin_temp.iloc[i, 0])
                        to_remove= [i for i in kin_temp.ID1 if i not in to_keep]
                        to_remove= list(set(to_remove))
                        remove= pd.DataFrame({'FID': to_remove})
                        remove['IID']= remove.FID
                        return remove
                df_list= list()
                flist= [infile for infile in input if 'pheno' in infile]
                for infile in flist:
                        x= pd.read_csv(infile, sep= '\t', header= 0)
                        df_list.append(x)
                d= pd.concat(df_list)
                pca= pd.read_csv(input[0], header= 0, sep= '\t')
                remove= selectUnrelated(d, d.IID)
                d= d[~d.IID.isin(remove)]
                d= pd.merge(d, pca, how= 'left', on= 'IID')
#                d['cohort']= d.cohort.astype('category').cat.codes
                d.to_csv(output[0], sep= '\t', header= True, index= False)



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

rule concat_heterozygosity:
	''
	input:
		expand('/mnt/work/pol/ROH/{cohort}/results/het/{{sample}}_excess_hom.txt', cohort= cohort_nms)
	output:
		'/mnt/work/pol/ROH/pheno/excess_{sample}.txt'
	run:
		df_list= list()
		for infile in input:
			d= pd.read_csv(infile, header= 0, sep= '\t')
			d['cohort']= np.where('harvestm12' in infile, 'harvestm12', np.where('harvestm24' in infile, 'harvestm24', np.where('rotterdam1' in infile, 'rotterdam1', np.where('rotterdam2' in infile, 'rotterdam2', np.where('normentfeb' in infile, 'normentfeb', 'normentmay')))))
			df_list.append(d)
		d= pd.concat(df_list)
		d.to_csv(output[0], header= True, index= False, sep= '\t')

rule replace_bp_cm_gap:
        'Replace UCSC gap bp position to cM.'
        input:
                '/mnt/work/pol/refdata/UCSC_gap.txt',
                '/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_combined_b37.txt',
		'/mnt/work/pol/ROH/1KG/1000GP_Phase3/chrX/genetic_map_chrX_nonPAR_combined_b37.txt'
        output:
                '/mnt/work/pol/ROH/1KG/cm_UCSC_gap.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                d.columns= ['chr', 'start', 'end', 'size', 'type']
		g= pd.read_csv(input[1], delim_whitespace= True, header= 0, names= ['chr', 'pos', 'rate', 'cM'])
                g= g[['chr', 'cM', 'pos']]
		g2= pd.read_csv(input[2], delim_whitespace= True, header= 0, names= ['pos', 'rate', 'cM'])
		g2['chr']= 23
		gm= pd.concat([g, g2])
                df_list= list()
                for CHR in set(d.chr):
                        temp_d= d.loc[d.chr== CHR, :]
                        temp_gm= gm.loc[gm.chr== CHR, :]
                        temp_d['cM1']= np.interp(temp_d.start, temp_gm['pos'], temp_gm['cM'])
                        temp_d['cM2']= np.interp(temp_d.end, temp_gm['pos'], temp_gm['cM'])
                        df_list.append(temp_d)
                x= pd.concat(df_list)
		x= x[['chr', 'cM1', 'cM2']]
                x.to_csv(output[0], sep= '\t', header= True, index= False)


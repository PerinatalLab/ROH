import pandas as pd
import numpy as np
import os
import gzip

cohort_nms= ['harvest','rotterdam1']
smpl_nms= ['maternal','paternal', 'fetal']
batch_nms= ['m12', 'm24']
CHR_nms= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

rule all:
	'Collect the main outputs of the workflow.'
	input:
#		expand('/mnt/work/pol/ROH/{cohort}/runs/frequency/ROH_frequency_{sample}', cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt', cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{cohort}/genotypes/maps/gene/{cohort}_{sample}_CHR{CHR}', cohort= cohort_nms, sample= smpl_nms, CHR= CHR_nms),
		expand('/mnt/work/pol/ROH/{cohort}/results/maps_cox/gene/cox_spont{sample}',cohort= cohort_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/{cohort}/pheno/IBD_parents.genome', cohort= cohort_nms)

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
                '~/soft/plink --bfile {params[0]} --indep-pairwise 50 5 0.9 --maf 0.05 --keep {input[1]} --make-bed --not-chr 23,24,25,26 --make-founders --out {params[1]}'

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
                '~/soft/plink --bfile {params[0]} --exclude range {input[2]} --indep-pairwise 50 5 0.9 --maf 0.05 --keep {input[1]} --make-bed --not-chr 23,24,25,26 --make-founders --out {params[1]}'

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

rule estimate_ROH:
        '''
        Obtain ROH estimates using PLINK 1.9.
        Configuration similar to :
                10.1186/1471-2164-12-460
        and used by:
                10.1371/journal.pgen.1007556
        '''
        input:
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.bed'
        output:
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom.indiv',
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom'
        params:
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}',
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}'
        shell:
                '''
                /home/pol.sole.navais/soft/plink --bfile {params[0]} --homozyg-window-snp 65 --homozyg-snp 65 --homozyg-kb 10 --homozyg-gap 500 --homozyg-window-missing 3 --homozyg-window-het 0 --homozyg-density 200 --out {params[1]}
                '''

rule rott_estimate_ROH:
        '''
        Obtain ROH estimates using PLINK 1.9.
        Configuration similar to :
                10.1186/1471-2164-12-460
        and used by:
                10.1371/journal.pgen.1007556
        '''
        input:
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.bed'
        output:
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom.indiv',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom'
        params:
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}'
        shell:
                '''
                /home/pol.sole.navais/soft/plink --bfile {params[0]} --homozyg-window-snp 65 --homozyg-snp 65 --homozyg-kb 10 --homozyg-gap 500 --homozyg-window-missing 3 --homozyg-window-het 0 --homozyg-density 200 --out {params[1]}
                '''

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
		'~/soft/plink --bfile {params[0]} --bmerge {params[1]} --out {params[2]}'

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
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm12_{sample}.bim',
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm24_{sample}.bim',
                '/mnt/work/pol/ROH/{cohort}/pheno/relatedness/relatedness_{sample}.kin0',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom',
                '/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom.indiv',
                '/mnt/work/pol/rotterdam1/pheno/rotterdam1_MFR.csv',
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.bim',
		'/mnt/work/pol/ROH/harvest/genotypes/prunedharvestm12_{sample}.fam'
        output:
                '/mnt/work/pol/ROH/{cohort}/pheno/runs_mfr_{sample}.txt'
        script:
                'scripts/pheno_file.py'

rule mapping_ROHs:
        'Obtain matrix (rows= position, columns = subject), with all ROHs per subject (1= homozygous part of ROH).'
        input:
                '/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom',
                '/mnt/work/pol/ROH/harvest/genotypes/prunedharvest{batch}_{sample}.bim',
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
                '/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_{sample}.bim',
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

rule parental_PLINK_rott:
	'Merge parental PLINK files and calculate IBD.'
	input:
		expand('/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_maternal.{ext}', ext= ['bed','bim','fam']),
		expand('/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_paternal.{ext}', ext= ['bed','bim','fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/rotterdam1/genotypes/parental/pruned_parental.{ext}', ext= ['bed','bim', 'fam', 'log']))
	params:
		'/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_maternal',
		'/mnt/work/pol/ROH/rotterdam1/genotypes/prunedrotterdam1_paternal',
		'/mnt/work/pol/ROH/rotterdam1/genotypes/parental/pruned_parental'
	shell:
		'~/soft/plink --bfile {params[0]} --bmerge {params[1]} --out {params[2]}'

rule parental_PLINK_harv:
	'Merge parental PLINK files and calculate IBD.'
	input:
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_maternal.{ext}', ext= ['bed','bim','fam']),
		expand('/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_paternal.{ext}', ext= ['bed','bim','fam'])
	output:
		temp(expand('/mnt/work/pol/ROH/harvest/genotypes/parental/pruned_parental.{ext}', ext= ['bed','bim', 'fam', 'log']))
	params:
		'/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_maternal',
		'/mnt/work/pol/ROH/harvest/genotypes/temp/prunedharvest_allbatch_paternal',
		'/mnt/work/pol/ROH/harvest/genotypes/parental/pruned_parental'
	shell:
		'~/soft/plink --bfile {params[0]} --bmerge {params[1]} --out {params[2]}'


rule mod_parental_fam_harvest:
	'Add family ID to parental PLINK fam file.'
	input:
		'/mnt/work/pol/ROH/harvest/genotypes/parental/pruned_parental.fam',
		'/mnt/work/pol/ROH/harvest/pheno/harvest_trios.txt'
	output:
		temp('/mnt/work/pol/ROH/harvest/genotypes/parental/mod_pruned_parental.fam')
	run:
		fam= pd.read_csv(input[0], sep= '\t', header= None)
		fam.columns= ['FID','IID','F', 'M', 'Sex','Pheno']
		d= pd.read_csv(input[1], sep= '\t')
		d.dropna(subset= ['Mother', 'Father'], inplace= True)
		d.drop_duplicates('Mother', inplace= True)
		d.drop_duplicates('Father', inplace= True)
		moms= d.loc[:, ['PREG_ID_1724','Mother']]
		dads= d.loc[:, ['PREG_ID_1724', 'Father']]
		dads.columns= ['FID2', 'Father']
		fam= pd.merge(fam, moms, right_on= 'Mother', left_on= 'IID', how= 'left')
		fam= pd.merge(fam, dads, right_on= 'Father', left_on= 'IID', how= 'left')
		fam['FID2']= np.where(~fam['Mother'].isnull(), fam.PREG_ID_1724, fam.FID2)
		fam['ng']= fam.groupby(fam.FID2.isnull()).cumcount() + 1 + max(fam.FID2)
		fam['FID2']= np.where(fam['FID2'].isnull(), fam.ng, fam.FID2)
		fam= fam[['FID2', 'IID', 'F', 'M', 'Sex', 'Pheno']]
		fam['FID2']= fam.FID2.astype(int)
		fam.to_csv(output[0], header=False, sep= '\t', index= False)

rule mod_parental_fam_rott:
	'Add family ID to parental PLINK fam file.'
	input:
		'/mnt/work/pol/ROH/rotterdam1/genotypes/parental/pruned_parental.fam',
		'/mnt/work/pol/ROH/rotterdam1/pheno/rotterdam1_trios.txt'
	output:
		temp('/mnt/work/pol/ROH/rotterdam1/genotypes/parental/mod_pruned_parental.fam')
	run:
		fam= pd.read_csv(input[0], sep= '\t', header= None)
		fam.columns= ['FID','IID','F', 'M', 'Sex','Pheno']
		d= pd.read_csv(input[1], sep= '\t')
		d.dropna(subset= ['Mother', 'Father'], inplace= True)
                d.drop_duplicates('Mother', inplace= True)
                d.drop_duplicates('Father', inplace= True)
		moms= d.loc[:, ['PREG_ID_315','Mother']]
		dads= d.loc[:, ['PREG_ID_315', 'Father']]
		dads.columns= ['FID2', 'Father']
		fam= pd.merge(fam, moms, right_on= 'Mother', left_on= 'IID', how= 'left')
		fam= pd.merge(fam, dads, right_on= 'Father', left_on= 'IID', how= 'left')
		fam['FID2']= np.where(~fam['Mother'].isnull(), fam.PREG_ID_315, fam.FID2)
		fam['ng']= fam.groupby(fam.FID2.isnull()).cumcount() + 1 + max(fam.FID2)
		fam['FID2']= np.where(fam['FID2'].isnull(), fam.ng, fam.FID2)
		fam= fam[['FID2', 'IID', 'F', 'M', 'Sex', 'Pheno']]
		fam['FID2']= fam.FID2.astype(int)
		fam.to_csv(output[0], header=False, sep= '\t', index= False)

rule plink_IBD:
	'Obtain IBD from PLINK parental duos.'
	input:
		'/mnt/work/pol/ROH/{cohort}/genotypes/parental/pruned_parental.bed',
		'/mnt/work/pol/ROH/{cohort}/genotypes/parental/pruned_parental.bim',
		'/mnt/work/pol/ROH/{cohort}/genotypes/parental/mod_pruned_parental.fam'
	output:
		'/mnt/work/pol/ROH/{cohort}/pheno/IBD_parents.genome'
	params:
		'/mnt/work/pol/ROH/{cohort}/pheno/IBD_parents'
	shell:
		'~/soft/plink --bed {input[0]} --bim {input[1]} --fam {input[2]} --genome rel-check --out {params[0]}'

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
		expand('/mnt/work/pol/ROH/1KG/1000GP_Phase3/genetic_map_chr{CHR}_combined_b37.txt', CHR= CHR_nms)
	shell:
		'''
		wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -P /mnt/work/pol/ROH/1KG/
		tar -xvzf /mnt/work/pol/ROH/1KG/1000GP_Phase3.tgz
		mv 1000GP_Phase3 /mnt/work/pol/ROH/1KG/
		rm /mnt/work/pol/ROH/1KG/1000GP_Phase3.tgz /mnt/work/pol/ROH/1KG/1000GP_Phase3/*hap.gz /mnt/work/pol/ROH/1KG/1000GP_Phase3/*.legend.gz /mnt/work/pol/ROH/1KG/1000GP_Phase3/1000GP_Phase3.sample
		'''

rule generate_report:
        'Generate report for harvest analysis.'
        input:
                expand('/mnt/work/pol/ROH/{{cohort}}/pheno/runs_mfr_{sample}.txt', sample= smpl_nms),
                '/mnt/work/pol/{cohort}/pheno/mod_{cohort}_q1_v9.csv',
                expand('/mnt/work/pol/ROH/harvest/runs/harvest{batch}_{sample}.hom', sample= smpl_nms, batch= batch_nms),
                expand('/mnt/work/pol/ROH/rotterdam1/runs/rotterdam1_{sample}.hom', sample= smpl_nms),
#                expand('/mnt/work/pol/ROH/rotterdam1/genotypes/maps/gene/rotterdam1_{sample}_CHR{CHR}', sample= smpl_nms, CHR= CHR_nms),
#                expand('/mnt/work/pol/ROH/harvest/genotypes/maps/gene/harvest{batch}_{sample}_CHR{CHR}', sample= smpl_nms, batch= batch_nms, CHR= CHR_nms),
#		expand('/mnt/work/pol/ROH/{{cohort}}/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}', CHR= CHR_nms, sample= smpl_nms)
        output:
                '/home/pol.sole.navais/ROH/reports/ROH_{cohort}_analysis.html'
        shell:
                """
                echo 'rmarkdown::render(input="scripts/ROH_harvest.Rmd", output_file="{output}")' | R --vanilla
                """



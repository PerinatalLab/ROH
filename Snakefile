import pandas as pd
import numpy as np
import os
import gzip

smpl_nms= ['maternal','paternal', 'fetal']
batch_nms= ['m12', 'm24']
CHR_nms= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

rule all:
	'Collect the main outputs of the workflow.'
	input:
		expand('/mnt/work/pol/ROH/runs/harvest_{batch}_{sample}.hom', batch= batch_nms, sample= smpl_nms),
		expand('/mnt/work/pol/ROH/pheno/runs_mfr_{sample}.txt', sample= smpl_nms),
		'/home/pol.sole.navais/ROH/reports/harvest_analysis.html',
		expand('/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz', sample= smpl_nms, CHR= CHR_nms),
		expand('/mnt/work/pol/ROH/runs/frequency/ROH_frequency_{sample}', sample= smpl_nms),
		expand('/mnt/work/pol/ROH/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}', sample= smpl_nms, CHR= CHR_nms)

rule ids_to_keep:
	'List of maternal, paternal and fetal ids acceptable by PLINK for --keep.'
	input:
		'/mnt/work/pol/harvest/pheno/harvest_linkage.csv'
	output:
		'/mnt/work/pol/ROH/pheno/maternal_ids',
		'/mnt/work/pol/ROH/pheno/paternal_ids',
		'/mnt/work/pol/ROH/pheno/fetal_ids'
	run:
		d= pd.read_csv(input[0], sep= ';')
		d['IID']= d.SentrixID_1
		d['FID']= d.IID
		d.drop(['PREG_ID_1724', 'FamilyID', 'SentrixID_2', 'SentrixID_3', 'Trio', 'noPhenoData', 'noGeneticData'], axis= 1, inplace= True)
		mat= d.loc[d.Role == 'Mother', ['FID', 'IID']]
		fet= d.loc[d.Role == 'Child', ['FID', 'IID']]
		fat= d.loc[d.Role == 'Father', ['FID', 'IID']]
		mat.to_csv(output[0], index= False, sep= '\t')
		fet.to_csv(output[2], index= False, sep= '\t')
		fat.to_csv(output[1], index= False, sep= '\t')

rule create_temp_folder:
	'Create a folder that will be temporary (removed when no longer in use)'
	output:
		temp(directory('/mnt/work/pol/ROH/genotypes/temp/'))
	shell:
		'mkdir {output}'

rule plink_split_bed:
        'Modify the bed file: remove CHR 23, 24, 25 and 26, maf <=0.01 and split file by sample.'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/{batch}/{batch}-genotyped.bed',
		'/mnt/work/pol/ROH/pheno/{sample}_ids'
        output:
                '/mnt/work/pol/ROH/genotypes/temp/harvest_{batch}_{sample}.bed',
                '/mnt/work/pol/ROH/genotypes/temp/harvest_{batch}_{sample}.prune.out'
        params:
                '/mnt/archive/HARVEST/delivery-fhi/data/genotyped/{batch}/{batch}-genotyped',
                '/mnt/work/pol/ROH/genotypes/temp/harvest_{batch}_{sample}'
        shell:
                '~/soft/plink --bfile {params[0]} --indep 50 5 10 --maf 0.01 --keep {input[1]} --make-bed --not-chr 23,24,25,26 --out {params[1]}'

rule plink_bfile_prune:
	'Exclude genetic variants in prune.out files (obtained with rule plink_split_bed).'
	input:
		'/mnt/work/pol/ROH/genotypes/temp/harvest_{batch}_{sample}.bed',
		'/mnt/work/pol/ROH/genotypes/temp/harvest_{batch}_{sample}.prune.out'
	output:
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}.bed',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}.bim',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}.fam'
	params:
		'/mnt/work/pol/ROH/genotypes/temp/harvest_{batch}_{sample}',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}'
	shell:
		'~/soft/plink --bfile {params[0]} --exclude {input[1]} --make-bed --out {params[1]}'

rule estimate_ROH:
        '''
        Obtain ROH estimates using PLINK 1.9.
        Configuration similar to :
                10.1186/1471-2164-12-460
        and used by:
                10.1371/journal.pgen.1007556
        '''
        input:
                '/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}.bed'
	output:
		'/mnt/work/pol/ROH/runs/harvest_{batch}_{sample}.hom.indiv',
		'/mnt/work/pol/ROH/runs/harvest_{batch}_{sample}.hom'
	params:
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}',
		'/mnt/work/pol/ROH/runs/harvest_{batch}_{sample}'
	shell:
		'''
		/home/pol.sole.navais/soft/plink --bfile {params[0]} --homozyg-window-snp 65 --homozyg-snp 65 --homozyg-kb 10 --homozyg-gap 500 --homozyg-window-missing 3 --homozyg-window-het 0 --homozyg-density 200 --out {params[1]} 
		'''
rule merge_m12_m24:
	'Merge PLINK files from the two batches for relatedness calculation.'
	input:
		'/mnt/work/pol/ROH/genotypes/prunedharvest_m12_{sample}.bed',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_m24_{sample}.bed'
	output:
		'/mnt/work/pol/ROH/genotypes/temp/prunedharvest_allbatch_{sample}.bed'
	params:
		'/mnt/work/pol/ROH/genotypes/prunedharvest_m12_{sample}',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_m24_{sample}',
		'/mnt/work/pol/ROH/genotypes/temp/prunedharvest_allbatch_{sample}'
	shell:
		'~/soft/plink --bfile {params[0]} --bmerge {params[1]} --out {params[2]}'

rule combine_pca:
	'Obtain pca for all samples.'
	input:
		'/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m12-founders-pca-covariates',
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m24-founders-pca-covariates',
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m12-offspring-pca-covariates',
                '/mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/m24-offspring-pca-covariates'
	output:
		'/mnt/work/pol/harvest/pheno/pca_fets_parents_m12_m24.txt'
	shell:
		'cat /mnt/archive/HARVEST/delivery-fhi/data/aux/pca-core/*-pca-covariates > {output}'			

rule relatedness:
	'Calculate relatedness using KING function from PLINK2.'
	input:
		'/mnt/work/pol/ROH/genotypes/temp/prunedharvest_allbatch_{sample}.bed',
		'/mnt/work/pol/ROH/pheno/{sample}_ids'
	output:
		'/mnt/work/pol/ROH/pheno/relatedness/relatedness_{sample}.kin0'
	params:
		'/mnt/work/pol/ROH/genotypes/temp/prunedharvest_allbatch_{sample}',
		'/mnt/work/pol/ROH/pheno/relatedness/relatedness_{sample}'
	shell:
		'~/soft/plink2 --bfile {params[0]} --keep {input[1]} --make-king-table --king-table-filter 0.03125 --out {params[1]}'

rule phenofile:
	'Merge all data necessary to create a phenotype file with ROH.'
	input:
		'/mnt/work/pol/ROH/runs/harvest_m12_{sample}.hom',
		'/mnt/work/pol/ROH/runs/harvest_m12_{sample}.hom.indiv',
		'/mnt/work/pol/ROH/runs/harvest_m24_{sample}.hom',
		'/mnt/work/pol/ROH/runs/harvest_m24_{sample}.hom.indiv',
		'/mnt/work/pol/harvest/pheno/harvest_mfr.csv',
		'/mnt/work/pol/harvest/pheno/harvest_linkage.csv',
		'/mnt/work/pol/harvest/pheno/pca_fets_parents_m12_m24.txt',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_m12_{sample}.bim',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_m24_{sample}.bim',
		'/mnt/work/pol/ROH/pheno/relatedness/relatedness_{sample}.kin0'
	output:
		'/mnt/work/pol/ROH/pheno/runs_mfr_{sample}.txt'
	script:
		'scripts/pheno_file.py'

rule make_report:
	'Generate report for harvest analysis.'
	input:
		expand('/mnt/work/pol/ROH/pheno/runs_mfr_{sample}.txt', sample= ['maternal', 'paternal', 'fetal'])
	output:
		'/home/pol.sole.navais/ROH/reports/harvest_analysis.html'
	shell:
		"""
		echo 'rmarkdown::render(input="scripts/ROH_harvest.Rmd", output_file="{output}")' | R --vanilla
		"""

rule mapping_ROHs:
	'Obtain matrix with all ROHs per subject with 1 if subject has ROH at specific position.'
	input:
		'/mnt/work/pol/ROH/runs/harvest_{batch}_{sample}.hom',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}.bim',
		'/mnt/work/pol/ROH/genotypes/prunedharvest_{batch}_{sample}.fam'
	output:
		'/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_{batch}_{sample}_chr{CHR}.txt'
	params:
		'/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_{batch}_{sample}_chr'
	threads: 6
	script:
		'scripts/map_ROHs.py'

rule merge_maps:
	'Merge maps from the two batches, one file per chromosome and sample.'
	input:
		'/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_m12_{sample}_chr{CHR}.txt',
		'/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_m24_{sample}_chr{CHR}.txt'
	output:
		'/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz'
	run:
		d12= pd.read_csv(input[0], sep= '\t', header= 0)
		d24= pd.read_csv(input[1], sep= '\t', header= 0)
		
		d= pd.merge(d12, d24, how= 'outer', on =['CHR', 'BP'])
		d.to_csv(output[0], sep= '\t', compression= 'gzip', index= False)
		os.remove(input[0])
		os.remove(input[1])

rule ROH_freq:
	'Count per-position relative frequency of ROHs.'
	input:
		expand('/mnt/work/pol/ROH/genotypes/maps/{{sample}}/maps_{{sample}}_chr{CHR}.txt.gz', CHR= CHR_nms)
	output:
		'/mnt/work/pol/ROH/runs/frequency/ROH_frequency_{sample}'
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
		'/mnt/work/pol/ROH/genotypes/maps/{sample}/maps_{sample}_chr{CHR}.txt.gz',
		'/mnt/work/pol/ROH/pheno/runs_mfr_{sample}.txt'
	output:
		'/mnt/work/pol/ROH/results/maps_cox/{sample}/cox_spont_{sample}_chr{CHR}'
	script:
		'scripts/cox_ROH.R'



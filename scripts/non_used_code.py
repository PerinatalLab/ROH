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



# Autozygosity and delivery timing in Norwegian family-trios  


**Parental genetic relatedness**, the sharing of two homologous alleles co-inherited from a common ancestor, increases the length of offspring’s genome covered by consecutive homozygous IBD (**autozygosity**). Long autozygous segments, a product of recent parental relatedness, are enriched in low-frequency and rare deleterious variants with recessive effects. If autozygosity is accurately detected, its mapping in genotyping array studies captures the effects of low-frequency and rare damaging non-genotyped variants that lie within these segments.  
While the genetic basis of most traits is largely shaped by an individual’s genome, the parental genome has a substantial role in shaping the genetic basis of other traits, such as **birth timing**. With a broad-sense heritability ranging between 25-40%, gestational duration maternal and dominance components account for 15 and 11% of its variance, respectively.  
Autozygosity can be detected through **runs-of-homozygosity** (ROH), but relies on the arbitrary selection of multiple parameters (minimum number of SNPs, distance between two consecutive ROHs, etc).  

Here, we leveraged the relation between parental relatedness and offspring ROH in 20000 family-trios, to systematically  select ROH calling parameters a posteriori and estimate the maternal, paternal and fetal effects of genome-wide autozygosity and to identify regions with recessive effects on spontaneous delivery risk.  


## Getting Started  

A single invocation of Snakemake, outputs the following reports:
	- Cohort-specific report for maternal, paternal and fetal effects of genome-wide and regional ROH on spontaneous delivery risk,  
	- Meta-analysis of different cohorts.

## Methods

In this study we used family-trio genetic data drawn from the Norwegian Mother and Child (MoBa), linked to the Medical Birth Registry of Norway (MFR). Only trios genotyped within the same batch were considered. Samples with major ethnicity other than CEU or cryptic relatedness to other samples were excluded. We also excluded pregnancies with a duration <154 (considered unviable) or ≥308 days, with missing gestational duration estimated by ultrasound, conceived by in-vitro fertilization, or affected by any of the following: abruptio placentae, placenta previa, polyhydramnion or congenital malformations.  

Autozygosity calling was performed using PLINK v.1.9, but calling parameters were selected a posteriori. In our framework, offspring autozygosity should be strongly associated with parental IBD. Steps followed for selecting ROH calling parameters:  
	1. Phase family-trio data using Eagle v.2. 
	2. Estimate the total length of shared IBD segments between parents (GERMLINE v.1.5.3).   
	3. Call ROHs in the offspring using different combinations of parameters:  
									- LD pruning threshold,  
									- genetic vs. physical distance,  
									- varying number of heterozyogtes allowed within ROHs, and  
									- minimum SNP count to call a ROH).  
	4. Calculate coefficient of determination between multiple offspring ROH and parental relatedness.
	5. Select the parameters that maximize the coefficient of determination between offspring ROH and parental relatedness.   
	6. Call ROHs in mothers, fathers and offspring – assuming no inbreeding differences between one generation.  





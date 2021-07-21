# About this repository 
The repository provides workflows described in the Common Workflow Language open standard, for "AFA": Computationally efficient Ancestral Frequency estimation in Admixed populations. 
"AFA" computes ancestry-specific allele frequencies in admixed populations using maximum likelihood estimation by modeling the conditional probability of having an allele given proportions of genetic ancestries. The method can use either local proportion ancestries, "LAFA", or global ancestry proportions "GAFA".  
# Contact information: 
Tamar Sofer:  tsofer@bwh.harvard.edu
Einat Granot-Hershkovitz:  egranot-hershkovitz@bwh.harvard.edu
# Required input for ancestry-global-proportion.cwl (GAFA):
-	plist: A data.frame with two columns (Split each gds file into blocks of 3,000 variants to calculate number of blocks in each chromosome):
1.	chromosome number
2.	block number

-	global_ancestry: a data.frame with a subject ID column, and an additional column for each of the ancestral proportions.
-	chromosome_gds: a SeqArray gds file. 
-	Annotation: a data.frame with a subject ID column, and a genetic consent column. 
-	myIndex – a number for running a specific job (one block). 
-	unrel – a data.frame with a subject ID column that includes unrelated individuals only. 
LAFA: This workflow calculates the ancestry-specific allele frequencies of bi-allelic genetic variants in admixed populations, based on local proportion ancestries. This workflow is specific for a three-way admixed population but can be adapted to any number of ancestries. 
# Required input for ancestry-local-proportion.cwl (LAFA):
-	plist: A data.frame with two columns (Split each gds file into blocks of 3,000 variants to calculate number of blocks in each chromosome):
1.	chromosome number
2.	block number
-	chromosome_gds: a SeqArray gds file. 
-	myIndex – a number for running a specific job (one block). 
-	unrel – a data.frame with a subject ID column that includes unrelated individuals only. 
-	African local ancestry file: a data.frame, each row represents one individual, each column represents one African local ancestry interval, taking a value of 0/1/2. This file only includes individuals with genetic consent. 
-	European local ancestry file: a data.frame, each row represents one individual, each column represents one European local ancestry interval, taking a value of 0/1/2. This file only includes individuals with genetic consent. 
-	Amerindian local ancestry file: a data.frame, each row represents one individual, each column represents one Amerindian local ancestry interval, taking a value of 0/1/2. This file only includes individuals with genetic consent. 
-	Annotation: an RData file with information on the Local ancestry intervals:
1.	chromosome: Chromosome no.
2.	pos.start 
3.	pos.end	
4.	snpID

-	chain file- liftover from hg38 to hg19 (in order to match each variant to its local ancestry interval, using the same assembly). 

#Output file column annotation
Each of the files has the following  columns: 
- CHR:	Chromosome number
- POS:	Chromosomal position (hg38)
- allele_a:	Reference allele
- allele_b:	Alternative allele
- Africa_estimated_freq:	Estimated African frequency 
- Africa_low_CI:	Low point of the 95% CI of the African frequency
- Africa_high_CI:	High point of the 95% CI of the African frequency
- Europe_estimated_freq:	Estimated European frequency
- Europe_low_CI:	Low point of the 95% CI of the European frequency 
- Europe_high_CI:	High point of the 95% CI of the European frequency
- America_estimated_freq:	Estimated Amerindian frequency 
- America_low_CI:	Low point of the 95% CI of the Amerindian frequency
- America_high_CI:	High point of the 95% CI of the Amerindian frequency
- boundary:	Boundary condition used when the function converges. 
- imputed:	Imputed genotype
- R2:	Squared correlation between input genotypes and imputed dosages


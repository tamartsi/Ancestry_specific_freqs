# About this repository 
A repository provide code for computing ancestry-specific allele frequencies in admixed populations using maximum likelihood estimation, and results from employing the code to estimates European, African, and Amerindian allele frequencies using Hispanic/Latino genomes from the Hispanic Community Health Study/Study of Latinos. 

# Code
Functions to compute ancestry-specific allele frequencies, as well as code for simulations, and some tests/examples, is provided in the Code folder. 

# HCHS_SOL_GAFA_estimates
This folder provides compressed files with ancestry-specific allele frequencies estimated using the "GAFA" method, which uses global ancestry estimates: estimated proportion of genomes inherited from each of the ancestral population used. The estimates are based on imputed data from the HCHS/SOL. See folder for information about columns. 


# HCHS_SOL_LAFA_estimates
This folder provides compressed files with ancestry-specific allele frequencies estimated using the "LAFA" method, which uses local ancestry estimates: for each of over 14,000 genomic intervals annotated, each person was assigned, using RFMix ancestral origin of these intervals in their two chromosome (or one chromosomal copy in the case of chromosome X). These are used as proportions of genetic ancestry in the intervals. The estimates are based on imputed data from the HCHS/SOL. See folder for information about columns. 

# Description of files providing ancestry-specific allele frequencies
Each files has a few columns, as follows: 
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

# Acknowledgements
We thank the staff and participants of the Hispanic Community Health Study/Study of Latinos for their important contributions.
